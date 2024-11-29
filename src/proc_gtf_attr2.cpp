#include <Rcpp.h>
#include <vector>
#include <string>
#include <regex>

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::string> processStrings(std::vector<std::string> x) {
    std::vector<std::string> result;
    std::regex quote("\"");
    std::regex key_value_split(" ");
    std::regex trailing_semicolon(";$");

    for (auto &str : x) {
        // Step 1: Remove the trailing semicolon
        str = std::regex_replace(str, trailing_semicolon, "");

        // Step 2: Split by "; "
        std::vector<std::string> attributes;
        size_t pos = 0;
        std::string delimiter = "; ";
        while ((pos = str.find(delimiter)) != std::string::npos) {
            attributes.push_back(str.substr(0, pos));
            str.erase(0, pos + delimiter.length());
        }
        attributes.push_back(str); // Add the last part

        // Step 3 and 4: Split by first space and remove quotes
        for (auto &attr : attributes) {
            size_t space_pos = attr.find(" ");
            if (space_pos != std::string::npos) {
                std::string key = attr.substr(0, space_pos);
                std::string value = attr.substr(space_pos + 1);
                value = std::regex_replace(value, quote, ""); // Remove quotes
                //result.push_back(key + " " + value);
                result.push_back(key);
                result.push_back(value);
            }
        }
    }
    return result;
}
