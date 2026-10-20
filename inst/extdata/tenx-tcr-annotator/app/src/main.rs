use amino::nucleotide_to_aminoacid_sequence;
use anyhow::{Context, Result, bail, ensure};
use debruijn::dna_string::DnaString;
use debruijn::Mer;
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use vdj_ann::annotate::{ContigAnnotation, Region};
use vdj_ann::refx::RefData;
use vdj_ann::vdj_features::{cdr1_start, cdr2_start, fr1_start, fr2_start, fr3_start};
use vdj_types::VdjRegion;

const HELP: &str = "tenx-tcr-annotate --input sequences.fasta --output annotations.csv\n\
    [--reference /path/to/regions.fa] [--json annotations.json]\n\
Uses Cell Ranger 10.1.0 annotation code on assembled DNA sequences.\n\
Coordinates are zero-based, end-exclusive, on the oriented output sequence.\n\
Only unambiguous A/C/G/T sequences are accepted. Cell/UMI/read metadata\n\
and clonotype inference are not produced by sequence-only annotation.\n";

fn parse_args() -> Result<BTreeMap<String, String>> {
    let mut parsed = BTreeMap::new();
    let mut args = std::env::args().skip(1);
    while let Some(key) = args.next() {
        if key == "--help" || key == "-h" {
            print!("{HELP}");
            std::process::exit(0);
        }
        ensure!(["--input", "--output", "--reference", "--json"].contains(&key.as_str()), "Unknown option: {key}\n{HELP}");
        let value = args.next().with_context(|| format!("Missing value for {key}"))?;
        ensure!(parsed.insert(key.clone(), value).is_none(), "Repeated option: {key}");
    }
    ensure!(parsed.contains_key("--input") && parsed.contains_key("--output"), "{HELP}");
    Ok(parsed)
}

fn read_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    let text = fs::read_to_string(path).with_context(|| format!("Read {}", path.display()))?;
    let mut records: Vec<(String, String)> = Vec::new();
    let mut ids = HashSet::new();
    for line in text.lines().map(str::trim).filter(|line| !line.is_empty()) {
        if let Some(header) = line.strip_prefix('>') {
            let id = header.split_whitespace().next().context("Empty FASTA header")?;
            ensure!(ids.insert(id.to_string()), "Duplicate FASTA identifier: {id}");
            records.push((id.to_string(), String::new()));
        } else {
            let (id, sequence) = records.last_mut().context("FASTA sequence before first header")?;
            ensure!(line.bytes().all(|b| b"ACGTacgt".contains(&b)), "Sequence {id} contains a base other than A/C/G/T; this wrapper rejects ambiguous bases rather than allowing the upstream packed-DNA representation to change them");
            sequence.push_str(&line.to_ascii_uppercase());
        }
    }
    ensure!(!records.is_empty(), "No FASTA records found");
    for (id, seq) in &records {
        ensure!(seq.len() >= 30, "Sequence {id} is shorter than 30 nucleotides");
    }
    Ok(records)
}

fn core_annotation(id: &str, dna: &DnaString, reference: &RefData) -> ContigAnnotation {
    // Sequence-only inputs have no quality, cell, read or UMI evidence. These
    // placeholders are required by the upstream API, and removed from our JSON.
    // The upstream API extracts a barcode before the first underscore and
    // panics if one is absent. FASTA identifiers need not follow that naming
    // convention; supply a temporary internal name and restore the original.
    let internal_id = if id.contains('_') { id.to_string() } else { format!("{id}_contig") };
    let mut annotation = ContigAnnotation::from_seq(dna, &vec![0; dna.len()], &internal_id, reference, 0, 0, false, false, None);
    annotation.contig_name = id.to_string();
    annotation
}

fn orientation_score(ann: &ContigAnnotation) -> (bool, usize, usize) {
    let vdj_count = ann.annotations.iter().filter(|a| matches!(a.feature.region_type, VdjRegion::V | VdjRegion::J)).count();
    let covered = ann.annotations.iter().map(|a| a.contig_match_end - a.contig_match_start).sum();
    (ann.is_productive(), vdj_count, covered)
}

fn region(sequence: &str, start: usize, stop: usize, origin: isize) -> Option<Region> {
    let frame = origin.rem_euclid(3) as usize;
    if start >= stop || stop > sequence.len() || start < frame {
        return None;
    }
    // Match Cell Ranger's region extraction: nucleotide endpoints may include
    // a final partial codon, whereas amino-acid endpoints are rounded down.
    let translated = nucleotide_to_aminoacid_sequence(sequence.as_bytes(), frame);
    let aa_start = (start - frame) / 3;
    let aa_stop = (stop - frame) / 3;
    if aa_stop > translated.len() {
        return None;
    }
    Some(Region {
        start,
        stop,
        nt_seq: sequence[start..stop].to_string(),
        aa_seq: String::from_utf8(translated[aa_start..aa_stop].to_vec()).ok()?,
    })
}

fn populate_regions(ann: &mut ContigAnnotation, reference: &RefData) -> Vec<String> {
    let mut notes = Vec::new();
    let Some(v) = ann.get_region(VdjRegion::V).cloned() else {
        notes.push("No V alignment; upstream framework/CDR regions unavailable".into());
        return notes;
    };
    let origin = v.contig_match_start as isize - v.annotation_match_start as isize;
    if let (Some(start), Some(j)) = (ann.cdr3_stop, ann.get_region(VdjRegion::J)) {
        ann.fwr4 = region(&ann.sequence, start, j.contig_match_end, origin);
    }
    if v.cigar.contains('I') || v.cigar.contains('D') {
        notes.push("V alignment contains an indel; FR1/CDR1/FR2/CDR2/FR3 transfer omitted".into());
        return notes;
    }
    let chain = v.feature.chain.to_string();
    if !["TRA", "TRB"].contains(&chain.as_str()) {
        notes.push("Framework/CDR1/CDR2 extraction in this wrapper is limited to TRA/TRB".into());
        return notes;
    }
    let Some(idx) = reference.id.iter().position(|&id| id as usize == v.feature.feature_id) else {
        notes.push("V reference identifier not found".into());
        return notes;
    };
    let aa = nucleotide_to_aminoacid_sequence(&reference.refs[idx].to_ascii_vec(), 0);
    // These are the unchanged upstream feature functions, applied to the V
    // reference exactly as in enclone_process/src/populate_features.rs.
    let offsets = [
        fr1_start(&aa, &chain),
        cdr1_start(&aa, &chain, false),
        fr2_start(&aa, &chain, false),
        cdr2_start(&aa, &chain, false),
        fr3_start(&aa, &chain, false),
    ];
    let Some(offsets) = offsets.into_iter().collect::<Option<Vec<_>>>() else {
        notes.push("Upstream feature search did not find every V-region boundary".into());
        return notes;
    };
    let Some(cdr3) = ann.cdr3_start else {
        notes.push("No CDR3; cannot delimit framework three".into());
        return notes;
    };
    let mut bounds: Vec<isize> = offsets.into_iter().map(|x| origin + 3 * x as isize).collect();
    bounds.push(cdr3 as isize);
    if bounds.windows(2).any(|w| w[0] >= w[1]) {
        notes.push("Inconsistent upstream feature order; V-region transfer omitted".into());
        return notes;
    }
    let extracted: Vec<Option<Region>> = bounds.windows(2).map(|w| {
        if w[0] < 0 || w[1] < 0 { None } else { region(&ann.sequence, w[0] as usize, w[1] as usize, origin) }
    }).collect();
    ann.fwr1 = extracted[0].clone();
    ann.cdr1 = extracted[1].clone();
    ann.fwr2 = extracted[2].clone();
    ann.cdr2 = extracted[3].clone();
    ann.fwr3 = extracted[4].clone();
    if extracted.iter().any(Option::is_none) {
        notes.push("Some V-region features are outside the covered sequence".into());
    }
    notes
}

fn put_region(row: &mut BTreeMap<String, String>, name: &str, region: Option<&Region>) {
    if let Some(r) = region {
        row.insert(name.into(), r.aa_seq.clone());
        row.insert(format!("{name}_nt"), r.nt_seq.clone());
        row.insert(format!("{name}_start"), r.start.to_string());
        row.insert(format!("{name}_end"), r.stop.to_string());
    }
}

fn csv_row(id: &str, ann: &ContigAnnotation, reversed: bool, notes: &[String]) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("sequence_id".into(), id.into());
    let consensus = id.replace("_consensus_", "_consensus");
    if let Some((clonotype, suffix)) = consensus.rsplit_once("_consensus") {
        if clonotype.starts_with("clonotype") && suffix.parse::<usize>().is_ok() {
            row.insert("consensus_id".into(), consensus.clone());
            row.insert("clonotype_id".into(), clonotype.into());
        }
    }
    row.insert("length".into(), ann.sequence.len().to_string());
    row.insert("sequence".into(), ann.sequence.clone());
    row.insert("reverse_complemented".into(), reversed.to_string());
    row.insert("annotation_notes".into(), notes.join("; "));
    row.insert("full_length".into(), ann.is_full_length().to_string());
    if let Some(p) = ann.productive { row.insert("productive".into(), p.to_string()); }
    if let Some(chain) = ann.chain_type() { row.insert("chain".into(), chain.to_string()); }
    for (region, prefix) in [(VdjRegion::V, "v"), (VdjRegion::D, "d"), (VdjRegion::J, "j"), (VdjRegion::C, "c")] {
        if let Some(a) = ann.get_region(region) {
            row.insert(format!("{prefix}_gene"), a.feature.gene_name.clone());
            if matches!(region, VdjRegion::V | VdjRegion::J) {
                row.insert(format!("{prefix}_start"), a.contig_match_start.to_string());
                row.insert(format!("{prefix}_end"), a.contig_match_end.to_string());
            }
            if region == VdjRegion::V { row.insert("v_end_ref".into(), a.annotation_match_end.to_string()); }
            if region == VdjRegion::J { row.insert("j_start_ref".into(), a.annotation_match_start.to_string()); }
        }
    }
    for (name, region) in [("fwr1", &ann.fwr1), ("cdr1", &ann.cdr1), ("fwr2", &ann.fwr2), ("cdr2", &ann.cdr2), ("fwr3", &ann.fwr3), ("fwr4", &ann.fwr4)] {
        put_region(&mut row, name, region.as_ref());
    }
    if let Some(aa) = &ann.cdr3 { row.insert("cdr3".into(), aa.clone()); }
    if let Some(nt) = &ann.cdr3_seq { row.insert("cdr3_nt".into(), nt.clone()); }
    if let Some(start) = ann.cdr3_start { row.insert("cdr3_start".into(), start.to_string()); }
    if let Some(stop) = ann.cdr3_stop { row.insert("cdr3_end".into(), stop.to_string()); }
    row
}

fn main() -> Result<()> {
    let args = parse_args()?;
    let reference_path = args.get("--reference").map(PathBuf::from).unwrap_or_else(|| {
        std::env::current_exe().unwrap().parent().unwrap().join("reference/regions.fa")
    });
    let reference = RefData::from_fasta(&reference_path).with_context(|| format!("Load reference {}", reference_path.display()))?;
    ensure!(!reference.refs.is_empty(), "Reference has no entries");
    let records = read_fasta(Path::new(&args["--input"]))?;
    let headers = "clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,fwr1,fwr1_nt,cdr1,cdr1_nt,fwr2,fwr2_nt,cdr2,cdr2_nt,fwr3,fwr3_nt,cdr3,cdr3_nt,fwr4,fwr4_nt,reads,umis,v_start,v_end,v_end_ref,j_start,j_start_ref,j_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,cdr3_start,cdr3_end,fwr4_start,fwr4_end,sequence_id,sequence,reverse_complemented,annotation_notes".split(',').collect::<Vec<_>>();
    let mut rows = Vec::new();
    let mut json = Vec::new();
    for (id, sequence) in &records {
        let dna = DnaString::from_dna_string(sequence);
        let forward = core_annotation(id, &dna, &reference);
        let reverse = core_annotation(id, &dna.rc(), &reference);
        let reversed = orientation_score(&reverse) > orientation_score(&forward);
        let mut ann = if reversed { reverse } else { forward };
        let notes = populate_regions(&mut ann, &reference);
        rows.push(csv_row(id, &ann, reversed, &notes));
        let mut item = serde_json::to_value(&ann)?;
        let object = item.as_object_mut().unwrap();
        for key in ["barcode", "quals", "read_count", "umi_count", "high_confidence", "is_cell", "filtered", "is_gex_cell", "is_asm_cell", "clonotype", "info", "junction_support", "sample", "fraction_of_reads_for_this_barcode_provided_as_input_to_assembly"] {
            object.remove(key);
        }
        object.insert("reverse_complemented".into(), reversed.into());
        object.insert("annotation_notes".into(), serde_json::to_value(&notes)?);
        json.push(item);
    }
    // Annotate all records before opening output, so malformed input cannot
    // produce an apparently successful partial result.
    if Path::new(&args["--input"]) == Path::new(&args["--output"]) { bail!("Input and output paths must differ"); }
    let mut writer = csv::Writer::from_path(&args["--output"])?;
    writer.write_record(&headers)?;
    for row in rows {
        writer.write_record(headers.iter().map(|key| row.get(*key).map_or("", String::as_str)))?;
    }
    writer.flush()?;
    if let Some(path) = args.get("--json") {
        fs::write(path, serde_json::to_string_pretty(&json)?)?;
    }
    eprintln!("Annotated {} sequences using Cell Ranger 10.1.0 code and {}", records.len(), reference_path.display());
    Ok(())
}
