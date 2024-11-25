use extendr_api::prelude::*;

/// Process a vector of strings, splitting by semicolon and cleaning the attributes.
///
/// # Arguments
/// * `x` - A vector of strings to process.
///
/// # Returns
/// A list of vectors, where each inner vector contains the processed attributes of a string.

/// @export
#[extendr]
pub fn process_attr_col_rust(x: Vec<String>) -> Vec<String> {
    // Create an output list
    let mut result = Vec::new();

    for line in x.iter() {
        let mut current_result = Vec::new();

        // Split the line by semicolons
        for segment in line.split(';') {
            // Trim whitespace and remove quotes
            let cleaned = segment
                .replace("\"", "") // Remove quotes
                .trim() // Remove leading/trailing whitespace
                .to_string();

            // Split by space and push each part
            for part in cleaned.split_whitespace() {
                current_result.push(part.to_string());
            }
        }

        // Append the processed vector to the result list
        // result.push(current_result); this creates a list of vectors Vec<Vec<String>>
        result.append(&mut current_result);
    }

    result
}

// Extendr module initialization
extendr_module! {
    mod igsc;
    fn process_attr_col_rust;
}
