#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cctype>

namespace fs = std::filesystem;

// =================================================================================
// Control Parameters Structure
// =================================================================================
struct ControlParams {
    int maxcyc = 50;
    std::string tmpdir = "./tmp";
    bool keeptmp = false;
    bool debug = false;
    
    void print_debug_info() const {
        if (debug) {
            std::cout << "\n=== Control Parameters (Debug) ===\n";
            std::cout << "maxcyc  = " << maxcyc << std::endl;
            std::cout << "tmpdir  = " << tmpdir << std::endl;
            std::cout << "keeptmp = " << (keeptmp ? "true" : "false") << std::endl;
            std::cout << "debug   = " << (debug ? "true" : "false") << std::endl;
            std::cout << "===================================\n" << std::endl;
        }
    }
};

// =================================================================================
// Helper Functions
// =================================================================================
bool execute_command(const std::string& cmd, bool debug = false) {
    if (debug) {
        std::cout << "Executing: " << cmd << std::endl;
    }
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Error: Command failed with exit code " << ret << std::endl;
        return false;
    }
    return true;
}
void replace_all(std::string& str, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}
std::string read_file(const fs::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << path << std::endl;
        return "";
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}
bool write_file(const fs::path& path, const std::string& content) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not write to file " << path << std::endl;
        return false;
    }
    file << content;
    return true;
}
struct XYZ {
    int atom_count = 0;
    std::string comment;
    std::vector<std::string> symbols;
    std::vector<double> coords;
    bool read(const fs::path& path) {
        std::ifstream file(path);
        if (!file) return false;
        file >> atom_count;
        std::getline(file, comment); 
        std::getline(file, comment); 
        symbols.resize(atom_count);
        coords.resize(atom_count * 3);
        for (int i = 0; i < atom_count; ++i) {
            file >> symbols[i] >> coords[i * 3] >> coords[i * 3 + 1] >> coords[i * 3 + 2];
        }
        return true;
    }
    std::string get_geom_str() const {
        std::stringstream ss;
        for (int i = 0; i < atom_count; ++i) {
            ss << symbols[i] << "    "
               << std::fixed << std::setprecision(10) << coords[i * 3] << "    "
               << coords[i * 3 + 1] << "    "
               << coords[i * 3 + 2] << (i == atom_count - 1 ? "" : "\n");
        }
        return ss.str();
    }
};

// =================================================================================
// Enhanced Parser Classes
// =================================================================================

class InputParser {
public:
    std::string prog;
    fs::path geom_file;
    std::string inp_tmplt1;
    std::string inp_tmplt2;
    std::map<std::string, std::string> grp_tmplt;
    ControlParams control;  // New control parameters

    bool parse(const fs::path& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open input file " << path << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);
            if (line.empty()) continue;

            if (line.rfind("%", 0) == 0) {
                std::stringstream line_ss(line.substr(1));
                std::string directive, value;
                line_ss >> directive;
                std::transform(directive.begin(), directive.end(), directive.begin(),
                               [](unsigned char c){ return std::tolower(c); });

                if (directive == "prog") {
                    line_ss >> prog;
                } else if (directive == "geom") {
                    line_ss >> geom_file;
                } else if (directive == "control") {
                    parse_control_section(file);
                } else if (directive == "inptmplt1" || directive == "inptmplt2" || directive == "grptmplt") {
                    std::stringstream content_ss;
                    std::string section_line;
                    while (std::getline(file, section_line)) {
                        std::string trimmed_line = section_line;
                        trimmed_line.erase(0, trimmed_line.find_first_not_of(" \t\n\r"));
                        if (trimmed_line == "end") {
                            break;
                        }
                        content_ss << section_line << '\n';
                    }

                    if (directive == "inptmplt1") {
                        inp_tmplt1 = content_ss.str();
                    } else if (directive == "inptmplt2") {
                        inp_tmplt2 = content_ss.str();
                    } else if (directive == "grptmplt") {
                        parse_grp_tmplt(content_ss.str());
                    }
                }
            }
        }
        return true;
    }

private:
    void parse_control_section(std::ifstream& file) {
        std::string line;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);
            
            if (line == "end") {
                break;
            }
            
            if (line.empty()) continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value = line.substr(eq_pos + 1);
                
                // Trim whitespace
                key.erase(0, key.find_first_not_of(" \t\n\r"));
                key.erase(key.find_last_not_of(" \t\n\r") + 1);
                value.erase(0, value.find_first_not_of(" \t\n\r"));
                value.erase(value.find_last_not_of(" \t\n\r") + 1);
                
                // Convert key to lowercase
                std::transform(key.begin(), key.end(), key.begin(),
                               [](unsigned char c){ return std::tolower(c); });
                
                // Parse different parameter types
                if (key == "maxcyc") {
                    try {
                        control.maxcyc = std::stoi(value);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: Invalid maxcyc value '" << value 
                                  << "', using default " << control.maxcyc << std::endl;
                    }
                } else if (key == "tmpdir") {
                    control.tmpdir = value;
                } else if (key == "keeptmp") {
                    std::transform(value.begin(), value.end(), value.begin(),
                                   [](unsigned char c){ return std::tolower(c); });
                    control.keeptmp = (value == "true" || value == "yes" || value == "1");
                } else if (key == "debug") {
                    std::transform(value.begin(), value.end(), value.begin(),
                                   [](unsigned char c){ return std::tolower(c); });
                    control.debug = (value == "true" || value == "yes" || value == "1");
                } else {
                    std::cerr << "Warning: Unknown control parameter '" << key << "'" << std::endl;
                }
            }
        }
    }
    
    void parse_grp_tmplt(const std::string& content) {
        std::stringstream ss(content);
        std::string line;
        if (std::getline(ss, line)) {
            std::smatch match;
            if (std::regex_search(line, match, std::regex(R"(state1\s*=\s*(\w+))"))) {
                grp_tmplt["state1"] = match[1];
            }
        }
        if (std::getline(ss, line)) {
            std::smatch match;
            if (std::regex_search(line, match, std::regex(R"(state2\s*=\s*(\w+))"))) {
                grp_tmplt["state2"] = match[1];
            }
        }
    }
};

class ConfigParser {
public:
    using ConfigMap = std::map<std::string, std::map<std::string, std::string>>;
    ConfigMap config;
    bool parse(const fs::path& path) {
        std::ifstream file(path);
        if (!file.is_open()) return false;
        std::string line, current_section, current_key, multiline_content;
        bool in_multiline = false;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);
            if (line.empty() || line[0] == '#') continue;
            if (in_multiline) {
                multiline_content += line + '\n';
                if (line.find('\'') != std::string::npos) {
                    multiline_content.pop_back();
                    multiline_content.erase(multiline_content.find_last_not_of("'") + 1);
                    config[current_section][current_key] = multiline_content;
                    in_multiline = false;
                }
                continue;
            }
            if (line[0] == '[' && line.back() == ']') {
                current_section = line.substr(1, line.length() - 2);
            } else {
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = line.substr(0, eq_pos);
                    key.erase(key.find_last_not_of(" \t") + 1);
                    std::string value = line.substr(eq_pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    if (value.front() == '\'') {
                        if (value.back() == '\'' && value.length() > 1) {
                            config[current_section][key] = value.substr(1, value.length() - 2);
                        } else {
                            in_multiline = true;
                            current_key = key;
                            multiline_content = value.substr(1) + '\n';
                        }
                    } else {
                        config[current_section][key] = value;
                    }
                }
            }
        }
        return true;
    }
};


// =================================================================================
// Enhanced Iterator Class with Control Parameters
// =================================================================================
class BaneMECPIterator {
public:
    BaneMECPIterator(InputParser input, fs::path executable_path) 
        : m_input(std::move(input)), m_exec_path(std::move(executable_path)) {
        // Use control parameters from input
        m_max_iter = m_input.control.maxcyc;
        m_tmpdir = m_input.control.tmpdir;
        m_keeptmp = m_input.control.keeptmp;
        m_debug = m_input.control.debug;
    }
    
    bool run() {
        // Print debug info if requested
        m_input.control.print_debug_info();
        
        if (!initialize()) return false;
        
        for (int i = 1; i <= m_max_iter; ++i) {
            std::cout << "\n===================================================\n"
                      << " MECP Iteration " << i << " (max: " << m_max_iter << ")"
                      << "\n===================================================\n";
            if (!run_iteration(i)) {
                std::cerr << "Iteration " << i << " failed." << std::endl;
                cleanup_if_needed();
                return false;
            }
            if (check_convergence()) {
                std::cout << "\n>>> MECP optimization converged in " << i << " steps. <<<" << std::endl;
                
                // Copy final result to original directory before cleanup
                if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
                    try {
                        fs::copy_file("final.xyz", m_original_dir / (m_base_name + "_final.xyz"), fs::copy_options::overwrite_existing);
                    } catch (const fs::filesystem_error& e) {
                        std::cerr << "Warning: Could not copy final result: " << e.what() << std::endl;
                    }
                } else {
                    fs::rename("final.xyz", m_base_name + "_final.xyz");
                }
                
                cleanup_if_needed();
                return true;
            }
            fs::rename("new.xyz", "astep" + std::to_string(i + 1) + ".xyz");
        }
        std::cerr << ">>> MECP optimization failed to converge in " << m_max_iter << " steps. <<<" << std::endl;
        cleanup_if_needed();
        return false;
    }
    
private:
    int m_max_iter;
    std::string m_tmpdir;
    bool m_keeptmp;
    bool m_debug;
    InputParser m_input;
    ConfigParser m_config;
    XYZ m_initial_geom;
    fs::path m_exec_path;
    fs::path m_original_dir;  // Store original working directory
    std::string m_base_name;
    bool m_is_external = false;
    
    bool initialize() {
        m_base_name = m_input.geom_file.stem().string();
        
        // Store original directory
        m_original_dir = fs::current_path();
        
        // Create and change to temporary directory if specified
        if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
            try {
                fs::create_directories(m_tmpdir);
                fs::current_path(m_tmpdir);
                if (m_debug) {
                    std::cout << "Debug: Created and changed to temporary directory: " << m_tmpdir << std::endl;
                    std::cout << "Debug: Original directory: " << m_original_dir << std::endl;
                }
                
                // Copy necessary files to tmp directory
                fs::path geom_src = m_original_dir / m_input.geom_file;
                fs::path banemecp_src = m_original_dir / "baneMECP.f";
                
                if (fs::exists(geom_src)) {
                    fs::copy_file(geom_src, m_input.geom_file.filename());
                    if (m_debug) {
                        std::cout << "Debug: Copied geometry file to tmp directory" << std::endl;
                    }
                } else {
                    std::cerr << "Error: Geometry file not found: " << geom_src << std::endl;
                    return false;
                }
                
                if (fs::exists(banemecp_src)) {
                    fs::copy_file(banemecp_src, "baneMECP.f");
                    if (m_debug) {
                        std::cout << "Debug: Copied baneMECP.f to tmp directory" << std::endl;
                    }
                } else {
                    std::cerr << "Error: baneMECP.f not found: " << banemecp_src << std::endl;
                    return false;
                }
                
            } catch (const fs::filesystem_error& e) {
                std::cerr << "Warning: Failed to setup tmpdir '" << m_tmpdir 
                          << "': " << e.what() << ". Using current directory." << std::endl;
                m_tmpdir = ".";
                fs::current_path(m_original_dir);
            }
        }
        
        if (!find_and_parse_config()) return false;
        if (!m_initial_geom.read(m_input.geom_file.filename())) {
            std::cerr << "Error: Could not read initial geometry from " << m_input.geom_file.filename() << std::endl;
            return false;
        }
        cleanup_files();
        if (!compile_mecp_solver(m_initial_geom.atom_count)) {
             std::cerr << "Error: Failed to compile MECP solver." << std::endl;
             return false;
        }
        fs::copy(m_input.geom_file.filename(), "astep1.xyz");
        append_to_trajectory("astep1.xyz");
        return true;
    }
    
    void cleanup_if_needed() {
        // Copy final results back to original directory if we were working in tmpdir
        if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
            try {
                // Copy final results back
                if (fs::exists("opt.trj")) {
                    fs::copy_file("opt.trj", m_original_dir / "opt.trj", fs::copy_options::overwrite_existing);
                }
                if (fs::exists(m_base_name + "_final.xyz")) {
                    fs::copy_file(m_base_name + "_final.xyz", m_original_dir / (m_base_name + "_final.xyz"), fs::copy_options::overwrite_existing);
                }
                
                // Change back to original directory
                fs::current_path(m_original_dir);
                if (m_debug) {
                    std::cout << "Debug: Returned to original directory: " << m_original_dir << std::endl;
                }
            } catch (const fs::filesystem_error& e) {
                std::cerr << "Warning: Error during cleanup: " << e.what() << std::endl;
            }
        }
        
        if (!m_keeptmp) {
            if (m_debug) {
                std::cout << "Debug: Cleaning up temporary files (keeptmp=false)" << std::endl;
            }
            
            // Remove temporary directory and all its contents
            if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
                try {
                    fs::remove_all(m_tmpdir);
                    if (m_debug) {
                        std::cout << "Debug: Removed temporary directory: " << m_tmpdir << std::endl;
                    }
                } catch (const fs::filesystem_error& e) {
                    if (m_debug) {
                        std::cout << "Debug: Could not remove tmpdir: " << e.what() << std::endl;
                    }
                }
            } else {
                // If working in current directory, clean up files normally
                cleanup_files();
            }
        } else {
            if (m_debug) {
                std::cout << "Debug: Keeping temporary files (keeptmp=true)" << std::endl;
                if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
                    std::cout << "Debug: Temporary files are in: " << fs::absolute(m_tmpdir) << std::endl;
                }
            }
        }
    }
    
    bool find_and_parse_config() {
        if (m_input.prog == "external") {
            m_is_external = true;
            if (m_debug) {
                std::cout << "Debug: Program is 'external'. Config file is optional." << std::endl;
            }
        }
        if (m_input.prog.empty()) {
            std::cerr << "Error: Program name is missing in the input file. Check the %Prog directive." << std::endl;
            return false;
        }
        std::string conf_name = m_input.prog + ".conf";
        fs::path search_paths[] = {
            m_original_dir / conf_name,
            m_exec_path / conf_name,
            fs::path(getenv("HOME")) / ".bane" / "mecp" / conf_name
        };
        for (const auto& path : search_paths) {
            if (fs::exists(path)) {
                if (m_debug) {
                    std::cout << "Debug: Using config file: " << path << std::endl;
                }
                // Copy config file to working directory if we're in tmpdir
                if (m_tmpdir != "." && m_tmpdir != "./" && !m_tmpdir.empty()) {
                    fs::copy_file(path, conf_name, fs::copy_options::overwrite_existing);
                    return m_config.parse(conf_name);
                } else {
                    return m_config.parse(path);
                }
            }
        }
        if (!m_is_external) {
             std::cerr << "Error: Config file '" << conf_name << "' not found." << std::endl;
             return false;
        }
        return true;
    }
    
    void cleanup_files() {
        if (m_debug) {
            std::cout << "Debug: Cleaning up temporary files..." << std::endl;
        }
        const std::vector<std::string> patterns = {
            "astep*", "MECP.state", "convg.tmp", "new.xyz", "final.xyz", 
            "MECP.x", "MECP_temp.f", "*.log", "*.gms", "*.gjf", "*.chk"
        };
        for (const auto& pattern : patterns) {
            std::string prefix = pattern.substr(0, pattern.find('*'));
            std::string suffix = pattern.substr(pattern.find('*') + 1);
            if (pattern.find('*') == std::string::npos) {
                prefix = pattern;
                suffix = "";
            }
            for (const auto& entry : fs::directory_iterator(".")) {
                std::string filename = entry.path().filename().string();
                bool match = false;
                if (pattern.find('*') == std::string::npos) {
                    if (filename == pattern) match = true;
                } else {
                    if (filename.rfind(prefix, 0) == 0 && (suffix.empty() || (filename.length() >= suffix.length() && filename.substr(filename.length() - suffix.length()) == suffix))) {
                       match = true;
                    }
                }
                if(match) {
                    if (m_debug) {
                        std::cout << "Debug: Removing file: " << filename << std::endl;
                    }
                    fs::remove(entry.path());
                }
            }
        }
        if (fs::exists("opt.trj")) fs::remove("opt.trj");
    }
    
    bool compile_mecp_solver(int natom) {
        if (m_debug) {
            std::cout << "Debug: Compiling baneMECP.f for " << natom << " atoms." << std::endl;
        }
        std::string fortran_code = read_file("baneMECP.f");
        if (fortran_code.empty()) {
            std::cerr << "Error: baneMECP.f not found in the current directory." << std::endl;
            return false;
        }
        replace_all(fortran_code, "parameter (Natom=3)", "parameter (Natom=" + std::to_string(natom) + ")");
        if (!write_file("MECP_temp.f", fortran_code)) return false;
        std::string compile_cmd = "gfortran -O3 -march=native -ffast-math -ffixed-line-length-none MECP_temp.f -o MECP.x";
        return execute_command(compile_cmd, m_debug);
    }
    
    bool run_iteration(int step) {
        std::string step_name = "astep" + std::to_string(step);
        std::string xyz_file = step_name + ".xyz";
        if (!prepare_qm_input(step, step_name, xyz_file)) return false;
        if (!run_qm_calculations(step_name)) return false;
        if (m_is_external && m_input.grp_tmplt.empty()) {
            if (m_debug) {
                std::cout << "Debug: Assuming external script generates .grad files." << std::endl;
            }
            if (!fs::exists(step_name + ".grad1") || !fs::exists(step_name + ".grad2")) {
                std::cerr << "Error: External script did not generate required .grad files." << std::endl;
                return false;
            }
        } else {
            if (!extract_results(step_name, xyz_file)) return false;
        }
        return execute_command("./MECP.x " + step_name, m_debug);
    }
    
    bool prepare_qm_input(int step, const std::string& step_name, const std::string& xyz_file) {
        XYZ current_geom;
        if (!current_geom.read(xyz_file)) {
            std::cerr << "Error: Failed to read geometry for step " << step << std::endl;
            return false;
        }
        std::string geom_str = current_geom.get_geom_str();
        auto process_template = [&](std::string tmpl) {
            replace_all(tmpl, "[geom]", geom_str);
            replace_all(tmpl, "[name]", step_name);
            replace_all(tmpl, "[xyzfile]", xyz_file);
            if (step == 1) {
                tmpl = std::regex_replace(tmpl, std::regex(R"(\[.*?\])"), "");
            } else {
                tmpl = std::regex_replace(tmpl, std::regex(R"(\[(.*?)\])"), "$1");
            }
            return tmpl;
        };
        std::string input1 = process_template(m_input.inp_tmplt1);
        std::string input2 = process_template(m_input.inp_tmplt2);
        if (m_is_external) {
            if (!write_file(step_name + "_state1.sh", input1)) return false;
            if (!write_file(step_name + "_state2.sh", input2)) return false;
        } else {
            std::string suffix = m_config.config["main"]["INPUT_SUFFIX"];
            if (!write_file(step_name + "_state1." + suffix, input1)) return false;
            if (!write_file(step_name + "_state2." + suffix, input2)) return false;
        }
        return true;
    }
    
    bool run_qm_calculations(const std::string& step_name) {
        if (m_is_external) {
            if (!execute_command("bash " + step_name + "_state1.sh", m_debug)) return false;
            if (!execute_command("bash " + step_name + "_state2.sh", m_debug)) return false;
        } else {
            std::string env = m_config.config["main"]["env"];
            std::string run_cmd_template = m_config.config["main"]["RUN_CMD"];
            std::string input_suffix = m_config.config["main"]["INPUT_SUFFIX"];
            std::string output_suffix = m_config.config["main"]["OUTPUT_SUFFIX"];
            auto run_single_state = [&](int state) {
                std::string actual_inp = step_name + "_state" + std::to_string(state) + "." + input_suffix;
                std::string actual_out = step_name + "_state" + std::to_string(state) + "." + output_suffix;
                std::string cmd = run_cmd_template;
                replace_all(cmd, "${ACTUAL_INP}", actual_inp);
                replace_all(cmd, "${ACTUAL_OUT}", actual_out);
                return execute_command(env + "\n" + cmd, m_debug);
            };
            if (!run_single_state(1)) return false;
            if (!run_single_state(2)) return false;
        }
        return true;
    }
    
    bool extract_results(const std::string& step_name, const std::string& xyz_file) {
         XYZ current_geom;
        if (!current_geom.read(xyz_file)) return false;
        std::string output_suffix = m_is_external ? "log" : m_config.config["main"]["OUTPUT_SUFFIX"];
        auto extract_single_state = [&](int state_num) {
            std::string state_type = m_input.grp_tmplt.at(state_num == 1 ? "state1" : "state2");
            const auto& rules = m_config.config.at(state_type);
            std::string output_file = step_name + "_state" + std::to_string(state_num) + "." + output_suffix;
            std::string content = read_file(output_file);
            if (content.empty()) return false;
            std::smatch e_match;
            std::regex e_regex(rules.at("E"));
            if (!std::regex_search(content, e_match, e_regex) || e_match.size() < 2) {
                std::cerr << "Error: Could not extract energy from " << output_file << std::endl;
                return false;
            }
            std::string energy = e_match[1].str();
            std::string gard_locate = rules.at("GARD.Loacte");
            size_t gard_pos = content.find(gard_locate);
            if (gard_pos == std::string::npos) {
                std::cerr << "Error: Could not locate gradient block in " << output_file << std::endl;
                return false;
            }
            std::stringstream content_stream(content.substr(gard_pos));
            std::string line;
            std::getline(content_stream, line);
            int n_skip = std::stoi(rules.at("GARD.NLineSkip"));
            for (int i = 0; i < n_skip; ++i) std::getline(content_stream, line);
            std::vector<int> cols;
            std::stringstream ss(rules.at("GARD.TargetColumns"));
            std::string col_str;
            while(std::getline(ss, col_str, ',')) cols.push_back(std::stoi(col_str));
            std::string end_by = rules.at("GARD.EndBy");
            std::stringstream grad_file_content;
            grad_file_content << current_geom.atom_count << "\n";
            grad_file_content << std::fixed << std::setprecision(10) << std::stod(energy) << "\n";
            for (int i = 0; i < current_geom.atom_count; ++i) {
                std::getline(content_stream, line);
                if ((!end_by.empty() && line.find(end_by) != std::string::npos) || (end_by.empty() && line.find_first_not_of(" \t\n\r") == std::string::npos)) {
                    std::cerr << "Error: Premature end of gradient block." << std::endl;
                    return false;
                }
                std::stringstream line_ss(line);
                std::vector<std::string> tokens;
                std::string token;
                while(line_ss >> token) tokens.push_back(token);
                grad_file_content << current_geom.symbols[i] << "    ";
                grad_file_content << std::stod(tokens.at(cols[0]-1)) << "    ";
                grad_file_content << std::stod(tokens.at(cols[1]-1)) << "    ";
                grad_file_content << std::stod(tokens.at(cols[2]-1)) << "\n";
            }
            return write_file(step_name + ".grad" + std::to_string(state_num), grad_file_content.str());
        };
        if (!extract_single_state(1)) return false;
        if (!extract_single_state(2)) return false;
        return true;
    }

    bool check_convergence() {
        if (!fs::exists("convg.tmp")) {
             std::cerr << "Warning: convg.tmp not found. Assuming not converged." << std::endl;
            return false;
        }
        std::string status = read_file("convg.tmp");
        
        // Trim whitespace from status string
        status.erase(0, status.find_first_not_of(" \t\n\r"));
        status.erase(status.find_last_not_of(" \t\n\r") + 1);

        if (status == "CONVERGED") {
            if (!fs::exists("final.xyz")) {
                std::cerr << "Error: Convergence indicated, but final.xyz was not created." << std::endl;
                return false; // Treat as error/not converged
            }
            append_to_trajectory("final.xyz");
            return true;
        } else {
            append_to_trajectory("new.xyz");
            return false;
        }
    }
    
    void append_to_trajectory(const fs::path& xyz_path) {
        if (!fs::exists(xyz_path)) {
            std::cerr << "Warning: Cannot append to trajectory, file not found: " << xyz_path << std::endl;
            return;
        }
        std::string content = read_file(xyz_path);
        std::ofstream traj_file("opt.trj", std::ios::app);
        if (traj_file.is_open()) {
            traj_file << content;
        }
    }
};

// =================================================================================
// Main Function (Enhanced with better error messages)
// =================================================================================
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <inputfile.inp>" << std::endl;
        std::cerr << "\nInput file should contain:" << std::endl;
        std::cerr << "  %Prog <program_name>" << std::endl;
        std::cerr << "  %Geom <geometry_file.xyz>" << std::endl;
        std::cerr << "  %Control" << std::endl;
        std::cerr << "    maxcyc=50" << std::endl;
        std::cerr << "    tmpdir=./tmp" << std::endl;
        std::cerr << "    keeptmp=false" << std::endl;
        std::cerr << "    debug=false" << std::endl;
        std::cerr << "  end" << std::endl;
        std::cerr << "  %InpTmplt1 ... end" << std::endl;
        std::cerr << "  %InpTmplt2 ... end" << std::endl;
        std::cerr << "  %GrpTmplt ... end" << std::endl;
        return 1;
    }
    
    fs::path input_file(argv[1]);
    if (!fs::exists(input_file)) {
        std::cerr << "Error: Input file not found: " << input_file << std::endl;
        return 1;
    }
    
    InputParser input_parser;
    if (!input_parser.parse(input_file)) {
        std::cerr << "Error: Failed to parse input file." << std::endl;
        return 1;
    }
    
    fs::path executable_path = fs::absolute(argv[0]).parent_path();
    BaneMECPIterator iterator(std::move(input_parser), std::move(executable_path));
    
    if (iterator.run()) {
        std::cout << "\nbanemecp finished successfully." << std::endl;
        return 0;
    } else {
        std::cerr << "\nbanemecp terminated with an error." << std::endl;
        return 1;
    }
}