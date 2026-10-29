# FASTQ compression utilities for Bash and zsh.
#
# Load these function definitions into the current shell:
#   source /Users/chris/Documents/R_packages/igsc/inst/extdata/compress_fastq/fastq_compression.bash
#
# Requirements:
#   relabel_fastq and compress_fastq_sequences: pigz and awk on PATH.
#   decode_fastq_triplets: Python 3 on PATH, using only its standard library.
#
# Compression inputs are gzip-compressed, unwrapped ASCII FASTQ records with
# four lines per read. The compression functions assume this layout; they do
# not fully validate FASTQ. Restoration takes a compressed sequence-only file.
# Compression outputs are gzip files regardless of the supplied filename suffix.
# Restored FASTQ output is plain text, or gzip-compressed when ending in .gz.
# Existing output files are overwritten. A failed run may leave partial output.
# Errors go to stderr. Return status is zero on success and nonzero on failure.
#
# The sequence-only formats retain bases and read order, but discard original
# headers, separator lines, and quality strings. Keep the read length separately
# to recover read boundaries; it is not stored in either sequence-only format.

# relabel_fastq -- optionally rewrite read names and qualities, then recompress.
#
# Usage:
#   relabel_fastq INPUT OUTPUT [--rename] [--quality] [--level 0-9]
#                 [--reads N] [--threads N]
#
# Arguments:
#   INPUT          Readable gzip-compressed FASTQ file.
#   OUTPUT         Destination gzip file; must differ from INPUT.
#
# Options:
#   --rename       Replace headers with @read1, @read2, ... in input order.
#                  Numbering starts at one on every function call. Default: off.
#   --quality      Replace each quality line with repeated I characters, using
#                  the first sequence length for every read. Default: off.
#   --level N      pigz compression level, an integer from 0 to 9. Default: 6.
#   --reads N      Keep at most the first N reads; N must be positive. Default:
#                  all reads. If fewer exist, keep all available reads.
#   --threads N    Pass a positive thread count to both pigz processes. Default:
#                  pigz chooses its thread count. The awk step remains serial.
#
# Behavior:
#   With no options, decompress and recompress without rewriting records.
#   Sequence lines and separator lines are preserved. Equal read lengths are
#   assumed when using --quality; this function does not check those lengths.
#   --reads may stop decompression early; the rest of the input is not checked.
#
# Example:
#   relabel_fastq input.fastq.gz renamed.fastq.gz --rename --quality --level 6 --threads 8 --reads 10000
#
# Returns:
#   Zero on success, nonzero on invalid arguments or a processing failure.
#   FASTQ data is written to OUTPUT, not stdout.
relabel_fastq() (
    set -o pipefail

    if [[ $# -lt 2 ]]; then
        echo "Usage: relabel_fastq INPUT OUTPUT [--rename] [--quality] [--level 0-9] [--reads N] [--threads N]" >&2
        return 1
    fi

    local input="$1" output="$2"
    local rename=0 replace_quality=0 level=6 max_reads=0
    local -a thread_args=()
    local pipeline_codes
    shift 2

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --rename) rename=1 ;;
            --quality) replace_quality=1 ;;
            --threads)
                if [[ $# -lt 2 || ! "$2" =~ ^[1-9][0-9]*$ ]]; then
                    echo "--threads requires a positive integer." >&2
                    return 1
                fi
                thread_args=(-p "$2")
                shift
                ;;
            --level)
                if [[ $# -lt 2 || ! "$2" =~ ^[0-9]$ ]]; then
                    echo "--level requires an integer from 0 to 9." >&2
                    return 1
                fi
                level="$2"
                shift
                ;;
            --reads)
                if [[ $# -lt 2 || ! "$2" =~ ^[1-9][0-9]*$ ]]; then
                    echo "--reads requires a positive integer." >&2
                    return 1
                fi
                max_reads="$2"
                shift
                ;;
            *)
                echo "Unknown argument: $1" >&2
                return 1
                ;;
        esac
        shift
    done

    if [[ ! -r "$input" || "$input" -ef "$output" ]]; then
        echo "Input must be readable and output must be a different file." >&2
        return 1
    fi

    pigz -dc "${thread_args[@]}" -- "$input" |
        awk -v rename="$rename" \
            -v replace_quality="$replace_quality" \
            -v max_reads="$max_reads" '
            max_reads && NR > 4 * max_reads { exit }
            NR % 4 == 1 && rename {
                print "@read" ++read_number
                next
            }
            NR == 2 && replace_quality {
                quality = $0
                gsub(/./, "I", quality)
            }
            NR % 4 == 0 && replace_quality {
                print quality
                next
            }
            { print }
        ' |
        pigz "-$level" "${thread_args[@]}" > "$output" || {
            # Capture exit codes immediately, using the current shell's variable.
            pipeline_codes="${PIPESTATUS[*]}${pipestatus[*]}"

            # pigz may report a broken pipe as SIGPIPE (141) or EPIPE (32).
            [[ "$max_reads" != 0 && (
                "$pipeline_codes" == "141 0 0" ||
                "$pipeline_codes" == "32 0 0"
            ) ]]
        }
)

## using triplets option made the file bigger and takes long compute!

# compress_fastq_sequences -- store only concatenated read sequences in gzip.
#
# Usage:
#   compress_fastq_sequences INPUT OUTPUT [--read-length N] [--level 0-9]
#                            [--threads N] [--reads N] [--triplets]
#
# Arguments:
#   INPUT           Readable gzip-compressed FASTQ file.
#   OUTPUT          Destination gzip file; must differ from INPUT.
#
# Options:
#   --read-length N Require every selected read to have exactly N bases, where
#                   N is positive. Default: 98. Reject length mismatches.
#   --level N       pigz compression level, an integer from 0 to 9. Default: 6.
#   --threads N     Pass a positive thread count to both pigz processes. Default:
#                   pigz chooses its thread count. The awk step remains serial.
#   --reads N       Keep at most the first N reads; N must be positive. Default:
#                   all reads. If fewer exist, keep all available reads.
#   --triplets      Apply the reversible custom triplet encoding before pigz.
#                   Default: off; bases are stored directly without line breaks.
#
# Output and validation:
#   Without --triplets, decompression produces one continuous sequence string.
#   With --triplets, decompression produces a custom encoded string, which must
#   be decoded with decode_fastq_triplets to recover the bases.
#   Original headers and qualities are discarded. Read order is retained.
#   A trailing carriage return is removed from each sequence line. Incomplete
#   selected records are rejected. --reads may stop before the input reaches EOF.
#   The expected read length is not saved in the output; record it separately.
#
# Triplet format:
#   Non-overlapping A/C/G/T triplets in lexical order map to these symbols:
#     ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
#   Examples: AAA -> A, AAC -> B, AAG -> C, AAT -> D, TTT -> /.
#   Grouping continues across read boundaries. A ~ prefix followed by three
#   literal bases preserves triplets containing N, lowercase, or other symbols.
#   A final ~ followed by one or two literal bases preserves the trailing bases.
#   This is a custom DNA encoding. It is not a standard Base64 file format.
#   Fewer precompression characters do not guarantee a smaller gzip file.
#
# Examples:
#   compress_fastq_sequences input.fastq.gz reads.L98.seq.gz --read-length 98 --threads 8
#   compress_fastq_sequences input.fastq.gz reads.L98.seq3.gz --read-length 98 --triplets --threads 8
#   compress_fastq_sequences input.fastq.gz first10000.seq.gz --reads 10000
#
# Recover sequences or FASTQ, using the original read length:
#   pigz -dc -- reads.L98.seq.gz | fold -w 98
#   decode_fastq_triplets reads.L98.seq.gz restored.fastq --read-length 98
#   decode_fastq_triplets reads.L98.seq3.gz restored.fastq --read-length 98 --triplets
#
# Returns:
#   Zero on success, nonzero on invalid arguments, unexpected sequence length,
#   incomplete selected records, or a processing failure. Data goes to OUTPUT.
compress_fastq_sequences() (
    set -o pipefail

    if [[ $# -lt 2 ]]; then
        echo "Usage: compress_fastq_sequences INPUT OUTPUT [--read-length N] [--level 0-9] [--threads N] [--reads N] [--triplets]" >&2
        return 1
    fi

    local input="$1" output="$2"
    local read_length=98 level=6 max_reads=0 pipeline_codes option
    local triplets=0
    local -a thread_args=()
    shift 2

    while [[ $# -gt 0 ]]; do
        option="$1"
        case "$option" in
            --triplets)
                triplets=1
                shift
                ;;
            --read-length|--threads|--reads)
                if [[ $# -lt 2 || ! "$2" =~ ^[1-9][0-9]*$ ]]; then
                    echo "$option requires a positive integer." >&2
                    return 1
                fi
                case "$option" in
                    --read-length) read_length="$2" ;;
                    --threads) thread_args=(-p "$2") ;;
                    --reads) max_reads="$2" ;;
                esac
                shift 2
                ;;
            --level)
                if [[ $# -lt 2 || ! "$2" =~ ^[0-9]$ ]]; then
                    echo "--level requires an integer from 0 to 9." >&2
                    return 1
                fi
                level="$2"
                shift 2
                ;;
            *)
                echo "Unknown argument: $option" >&2
                return 1
                ;;
        esac
    done

    if [[ ! -r "$input" || "$input" -ef "$output" ]]; then
        echo "Input must be readable and output must be a different file." >&2
        return 1
    fi

    pigz -dc "${thread_args[@]}" -- "$input" |
        awk -v read_length="$read_length" -v max_reads="$max_reads" \
            -v triplets="$triplets" '
            BEGIN {
                # ACGT triplets in lexical order map to these 64 symbols.
                alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
                bases = "ACGT"
                for (a = 1; a <= 4; a++)
                    for (b = 1; b <= 4; b++)
                        for (c = 1; c <= 4; c++) {
                            triplet = substr(bases, a, 1) substr(bases, b, 1) substr(bases, c, 1)
                            code[triplet] = substr(alphabet, ++index_number, 1)
                        }
            }
            # Emit a read as raw bases or complete triplet codes. Global carry
            # retains up to two bases for the next read; END flushes the tail.
            # text, i, and triplet are local scratch parameters.
            function emit_sequence(sequence, text, i, triplet) {
                if (!triplets) {
                    printf "%s", sequence
                    return
                }
                text = carry sequence
                for (i = 1; i + 2 <= length(text); i += 3) {
                    triplet = substr(text, i, 3)
                    if (triplet in code) printf "%s", code[triplet]
                    else printf "~%s", triplet
                }
                carry = substr(text, i)
            }
            NR % 4 == 2 {
                sub(/\r$/, "")
                if (length($0) != read_length) {
                    print "Unexpected read length at line " NR > "/dev/stderr"
                    failed = 1
                    exit 1
                }
                sequence = $0
            }
            NR % 4 == 0 {
                emit_sequence(sequence)
                reads++
                if (max_reads && reads >= max_reads) exit
            }
            END {
                if (!failed && NR % 4 != 0) {
                    print "Incomplete FASTQ record" > "/dev/stderr"
                    exit 1
                }
                # Escape a final one or two bases instead of discarding them.
                if (!failed && triplets && length(carry)) printf "~%s", carry
            }
        ' |
        pigz "-$level" "${thread_args[@]}" > "$output" || {
            pipeline_codes="${PIPESTATUS[*]}${pipestatus[*]}"

            # Allow a broken pipe only when the read limit stops processing early.
            [[ "$max_reads" != 0 && (
                "$pipeline_codes" == "141 0 0" ||
                "$pipeline_codes" == "32 0 0"
            ) ]]
        }
)

# decode_fastq_triplets -- restore FASTQ from compressed sequence-only data.
#
# Usage:
#   decode_fastq_triplets INPUT OUTPUT [--read-length N] [--triplets|--no-triplets]
#
# Arguments:
#   INPUT          A gzip file created by compress_fastq_sequences, with or
#                  without --triplets. Original FASTQ input is not supported.
#   OUTPUT         Destination FASTQ file; must differ from INPUT. A filename
#                  ending in .gz selects gzip compression at level 6. Other
#                  filenames receive plain FASTQ. Existing files are overwritten.
#
# Options:
#   --read-length N Split the recovered sequence into reads of exactly N bases;
#                   N must be positive. Default: 98. Use the compression length.
#   --triplets      Decode the custom triplet format before splitting reads.
#                   Supply this if compression used --triplets. Default: off.
#   --no-triplets   Explicitly treat the decompressed input as raw bases. This
#                   is the default. If both flags appear, the last one wins.
#
# Output and validation:
#   Write four lines per read: @read1, @read2, ... as the header; the recovered
#   sequence; a literal + line; and I repeated to match the specified read length.
#   Names restart at one per call. Original names and qualities are not retained.
#   Encoding mode and read length are not stored in the input; set them to match
#   compression. There is no automatic encoding detection.
#   Processing uses bounded chunks and preserves sequence order and case.
#   Reject invalid triplet symbols or escapes, line breaks in sequence data,
#   and a final incomplete read. Empty input produces an empty FASTQ file.
#   File data goes to OUTPUT, not stdout. Failures may leave partial output.
#
# Requirements:
#   Python 3 on PATH, using its standard gzip, itertools, and sys modules.
#   gzip decompression, optional output compression, and formatting are serial.
#
# Examples:
#   decode_fastq_triplets reads.L98.seq.gz restored.fastq --read-length 98
#   decode_fastq_triplets reads.L98.seq3.gz restored.fastq --read-length 98 --triplets
#   decode_fastq_triplets reads.L98.seq3.gz restored.fastq.gz --triplets
#
# Returns:
#   Zero on success; nonzero for invalid arguments, unreadable or invalid gzip
#   input, invalid sequence data or triplet encoding, or output failures.
#   Errors go to stderr. Already written output is not rolled back on failure.
decode_fastq_triplets() (
    if [[ $# -lt 2 ]]; then
        echo "Usage: decode_fastq_triplets INPUT OUTPUT [--read-length N] [--triplets|--no-triplets]" >&2
        return 1
    fi

    local input="$1" output="$2" read_length=98 triplets=0
    shift 2

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --triplets) triplets=1; shift ;;
            --no-triplets) triplets=0; shift ;;
            --read-length)
                if [[ $# -lt 2 || ! "$2" =~ ^[1-9][0-9]*$ ]]; then
                    echo "--read-length requires a positive integer." >&2
                    return 1
                fi
                read_length="$2"
                shift 2
                ;;
            *)
                echo "Unknown argument: $1" >&2
                return 1
                ;;
        esac
    done

    if [[ ! -r "$input" || "$input" -ef "$output" ]]; then
        echo "Input must be readable and output must be a different file." >&2
        return 1
    fi

    python3 - "$input" "$output" "$read_length" "$triplets" <<'PY'
import gzip
import itertools
import sys


def sequence_chunks(source, triplets):
    """Yield raw bases, optionally decoding triplets and literal escapes."""
    alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
    table = dict(zip(alphabet, (bytes(t) for t in itertools.product(b"ACGT", repeat=3))))
    pending = b""
    while True:
        block = source.read(65536)
        if not block:
            break
        if not triplets:
            yield block
            continue
        data = pending + block
        decoded = bytearray()
        position = 0
        while position < len(data):
            symbol = data[position]
            if symbol == ord("~"):
                if len(data) - position < 4:
                    break
                decoded.extend(data[position + 1:position + 4])
                position += 4
            else:
                if symbol not in table:
                    raise ValueError(f"invalid encoded symbol: {bytes([symbol])!r}")
                decoded.extend(table[symbol])
                position += 1
        pending = data[position:]
        yield decoded
    if pending:
        if pending[0] != ord("~") or len(pending) not in (2, 3):
            raise ValueError("incomplete literal escape")
        yield pending[1:]


def write_fastq(chunks, output, read_length):
    """Write complete reads with sequential names and constant I qualities."""
    pending = b""
    quality = b"I" * read_length
    read_number = 0
    for chunk in chunks:
        data = pending + chunk
        if b"\n" in data or b"\r" in data:
            raise ValueError("sequence data contains line breaks")
        complete = len(data) // read_length * read_length
        records = bytearray()
        for position in range(0, complete, read_length):
            read_number += 1
            sequence = data[position:position + read_length]
            records.extend(b"@read" + str(read_number).encode("ascii") + b"\n")
            records.extend(sequence + b"\n+\n" + quality + b"\n")
        output.write(records)
        pending = data[complete:]
    if pending:
        raise ValueError(
            f"incomplete final read: {len(pending)} bases; expected {read_length}"
        )


try:
    input_path, output_path = sys.argv[1:3]
    read_length = int(sys.argv[3])
    triplets = sys.argv[4] == "1"
    with gzip.open(input_path, "rb") as source:
        if output_path.endswith(".gz"):
            destination = gzip.open(output_path, "wb", compresslevel=6)
        else:
            destination = open(output_path, "wb")
        with destination as output:
            write_fastq(sequence_chunks(source, triplets), output, read_length)
except (OSError, EOFError, ValueError) as error:
    print(f"decode_fastq_triplets: {error}", file=sys.stderr)
    sys.exit(1)
PY
)



# recover sequences: pigz -dc -- reads.L98.seq.gz | fold -w 98

# recover FASTQ without triplet decoding:
# decode_fastq_triplets reads.L98.seq.gz restored.fastq --read-length 98

# recover FASTQ with triplet decoding:
# decode_fastq_triplets reads.L98.seq3.gz restored.fastq --read-length 98 --triplets
