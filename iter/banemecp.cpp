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
    std::string tmpdir = ".";
    bool keeptmp = false;
    bool debug = false;
    bool restart = false;
    
    // 新增收敛判据参数
    double tde = 5e-5;        // 能量gap阈值
    double tgmax = 7e-4;      // 最大梯度阈值  
    double tgrms = 5e-4;      // RMS梯度阈值
    double tdxmax = 4e-3;     // 最大位移阈值
    double tdxrms = 2.5e-3;   // RMS位移阈值
    double stpmx = 0.1;       // 置信半径/最大步长
    
    void print_debug_info() const {
        if (debug) {
            std::cout << "\n=== Control Parameters (Debug) ===\n";
            std::cout << "maxcyc  = " << maxcyc << std::endl;
            std::cout << "tmpdir  = " << tmpdir << std::endl;
            std::cout << "keeptmp = " << (keeptmp ? "true" : "false") << std::endl;
            std::cout << "debug   = " << (debug ? "true" : "false") << std::endl;
            std::cout << "restart = " << (restart ? "true" : "false") << std::endl;
            std::cout << "=== Convergence Criteria ===\n";
            std::cout << "tde     = " << tde << std::endl;
            std::cout << "tgmax   = " << tgmax << std::endl;
            std::cout << "tgrms   = " << tgrms << std::endl;
            std::cout << "tdxmax  = " << tdxmax << std::endl;
            std::cout << "tdxrms  = " << tdxrms << std::endl;
            std::cout << "stpmx   = " << stpmx << std::endl;
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
        std::map<std::string, std::function<void(const std::string&)>> parsers = {
            {"maxcyc", [this](const std::string& v) { control.maxcyc = std::stoi(v); }},
            {"tmpdir", [this](const std::string& v) { control.tmpdir = v; }},
            {"keeptmp", [this](const std::string& v) { control.keeptmp = (v == "true"); }},
            {"debug", [this](const std::string& v) { control.debug = (v == "true"); }},
            {"restart", [this](const std::string& v) { control.restart = (v == "true"); }},
            {"tde", [this](const std::string& v) { control.tde = std::stod(v); }},
            {"tgmax", [this](const std::string& v) { control.tgmax = std::stod(v); }},
            {"tgrms", [this](const std::string& v) { control.tgrms = std::stod(v); }},
            {"tdxmax", [this](const std::string& v) { control.tdxmax = std::stod(v); }},
            {"tdxrms", [this](const std::string& v) { control.tdxrms = std::stod(v); }},
            {"stpmx", [this](const std::string& v) { control.stpmx = std::stod(v); }}
        };

        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);

            if (line == "end") break;
            if (line.empty()) continue;

            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value = line.substr(eq_pos + 1);

                // Trim
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                std::transform(key.begin(), key.end(), key.begin(), ::tolower);
                // 注意：value不要转换为小写，因为数值不需要

                if (parsers.count(key)) {
                    try {
                        parsers[key](value);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: Failed to parse " << key << " = " << value 
                                  << ", using default value. Error: " << e.what() << std::endl;
                    }
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
        m_restart = m_input.control.restart;
    }
    
    bool run() {
        // Print debug info if requested
        m_input.control.print_debug_info();
        
        if (!initialize()) return false;
        
        for (int i = m_current_step; i <= m_max_iter; ++i) {
            std::cout << "\n===================================================\n"
                      << " MECP Iteration " << i << " (max: " << m_max_iter << ")"
                      << "\n===================================================\n";
            std::cout.flush();
            if (!run_iteration(i)) {
                std::cerr << "Iteration " << i << " failed." << std::endl;
                cleanup_if_needed();
                return false;
            }
            if (check_convergence()) {
                std::cout << "\n>>> MECP optimization converged in " << i << " steps. <<<" << std::endl;
                
                // Copy final result to original directory before cleanup
                if (using_tmpdir()) {
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
    bool m_restart;       // New flag: restart control
    int m_current_step = 1; // New flag: current step number
    InputParser m_input;
    ConfigParser m_config;
    XYZ m_initial_geom;
    fs::path m_exec_path;
    fs::path m_original_dir;  // Store original working directory
    std::string m_base_name;
    bool m_is_external = false;
    
    // New function: check if InpTmplt2 is provided
    bool has_inp_tmplt2() const {
        return !m_input.inp_tmplt2.empty();
    }
    
    // New function: read step number from convg.tmp
    int read_step_from_convg() {
        if (!fs::exists("convg.tmp")) {
            if (m_debug) {
                std::cout << "Debug: convg.tmp not found for restart step detection" << std::endl;
            }
            return -1;
        }
        
        std::string content = read_file("convg.tmp");
        if (content.empty()) {
            if (m_debug) {
                std::cout << "Debug: convg.tmp is empty" << std::endl;
            }
            return -1;
        }
        
        // Parse step number, e.g. Step: 37
        std::regex step_regex(R"(Step:\s*(\d+))");
        std::smatch match;
        
        if (std::regex_search(content, match, step_regex)) {
            int step = std::stoi(match[1].str());
            if (m_debug) {
                std::cout << "Debug: Found step " << step << " in convg.tmp" << std::endl;
            }
            return step;
        }
        
        if (m_debug) {
            std::cout << "Debug: Could not parse step number from convg.tmp" << std::endl;
        }
        return -1;
    }
    
    // New function: selective file cleanup
    void cleanup_files_selective(bool preserve_restart_files = false) {
        if (m_debug) {
            std::cout << "Debug: Selective cleanup (preserve_restart_files=" 
                      << (preserve_restart_files ? "true" : "false") << ")" << std::endl;
        }
        
        std::vector<std::string> patterns_to_clean = {
            "MECP.x", "MECP_temp.f", "*.log", "*.gms", "*.gjf", "*.chk", "new.xyz", "final.xyz"
        };
        
        if (!preserve_restart_files) {
            // If not preserving restart files, then clean all temp files
            patterns_to_clean.insert(patterns_to_clean.begin(), 
                {"astep*", "MECP.state", "convg.tmp", "opt.trj"});
        }
        
        for (const auto& pattern : patterns_to_clean) {
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
    }
    
    bool using_tmpdir() const {
        return !m_tmpdir.empty() && m_tmpdir != ".";
    }

    bool initialize() {
        m_base_name = m_input.geom_file.stem().string();
        
        // Print information about InpTmplt2
        if (m_debug) {
            std::cout << "Debug: InpTmplt2 " << (has_inp_tmplt2() ? "provided" : "not provided") << std::endl;
            if (!has_inp_tmplt2()) {
                std::cout << "Debug: Second state calculation will be skipped" << std::endl;
            }
        }
        
        // Store original directory
        m_original_dir = fs::current_path();

        // Create and change to temporary directory if specified
        if (using_tmpdir()) {
            try {
                // Remove existing tmpdir if it exists to avoid conflicts
                if (fs::exists(m_tmpdir)) {
                    if (m_debug) {
                        std::cout << "Debug: Removing existing tmpdir: " << m_tmpdir << std::endl;
                    }
                    fs::remove_all(m_tmpdir);
                }

                // Create fresh tmpdir and change to it
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
                    fs::copy_file(geom_src, m_input.geom_file.filename(), fs::copy_options::overwrite_existing);
                    if (m_debug) {
                        std::cout << "Debug: Copied geometry file to tmp directory" << std::endl;
                    }
                } else {
                    std::cerr << "Error: Geometry file not found: " << geom_src << std::endl;
                    return false;
                }

                if (fs::exists(banemecp_src)) {
                    fs::copy_file(banemecp_src, "baneMECP.f", fs::copy_options::overwrite_existing);
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
        
        // Check if we should restart based on existing files and restart flag
        bool has_state_file = fs::exists("MECP.state");
        bool has_convg_file = fs::exists("convg.tmp");
        bool should_restart = m_restart && (has_state_file || has_convg_file);

        if (should_restart) {
            int last_step = read_step_from_convg();
            if (last_step > 0) {
                m_current_step = last_step + 1;
                std::cout << ">>> Restarting from step " << m_current_step << " <<<" << std::endl;
            } else {
                m_current_step = 1;
            }
            // Selective cleanup, preserve restart-related files
            cleanup_files_selective(true);
        } else {
            if (has_state_file && !m_restart) {
                std::cout << ">>> Found MECP.state but restart=false. Starting fresh... <<<" << std::endl;
            } else if (!has_state_file && !has_convg_file) {
                std::cout << ">>> No restart files found. Starting fresh optimization... <<<" << std::endl;
            }
            
            // Full cleanup
            cleanup_files_selective(false);
            m_current_step = 1;
        }
        
        if (!compile_mecp_solver(m_initial_geom.atom_count)) {
             std::cerr << "Error: Failed to compile MECP solver." << std::endl;
             return false;
        }
        
        // If not restart or first step, prepare first step
        if (!should_restart || m_current_step == 1) {
            fs::copy(m_input.geom_file.filename(), "astep1.xyz");
            append_to_trajectory("astep1.xyz");
        }
        
        return true;
    }
    
    void cleanup_if_needed() {
        if (using_tmpdir()) {
            try {
                // Copy results back
                std::vector<std::string> files_to_copy = {"opt.trj", m_base_name + "_final.xyz"};
                for (const auto& file : files_to_copy) {
                    if (fs::exists(file)) {
                        fs::copy_file(file, m_original_dir / file, fs::copy_options::overwrite_existing);
                    }
                }
                
                fs::current_path(m_original_dir);
                
                if (!m_keeptmp) {
                    fs::remove_all(m_tmpdir);
                }
            } catch (const fs::filesystem_error& e) {
                std::cerr << "Warning: Cleanup error: " << e.what() << std::endl;
            }
        } else if (!m_keeptmp) {
            cleanup_files_selective(false);
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
                bool parse_success;
                if (using_tmpdir()) {
                    fs::copy_file(path, conf_name, fs::copy_options::overwrite_existing);
                    parse_success = m_config.parse(conf_name);
                } else {
                    parse_success = m_config.parse(path);
                }
                
                if (parse_success) {
                    // Check if RUN_CMD exists in [main] section
                    if (m_config.config.count("main") == 0 || m_config.config["main"].count("RUN_CMD") == 0) {
                        if (m_debug) {
                            std::cout << "Debug: RUN_CMD not found in config. Switching to external mode." << std::endl;
                        }
                        m_is_external = true;
                    }
                }
                return parse_success;
            }
        }
        if (!m_is_external) {
             std::cerr << "Error: Config file '" << conf_name << "' not found." << std::endl;
             return false;
        }
        return true;
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
            if (!fs::exists(step_name + ".grad1")) {
                std::cerr << "Error: External script did not generate required .grad1 file." << std::endl;
                return false;
            }
            // Only check for .grad2 if InpTmplt2 is provided
            if (has_inp_tmplt2() && !fs::exists(step_name + ".grad2")) {
                std::cerr << "Error: External script did not generate required .grad2 file." << std::endl;
                return false;
            }
        } else {
            if (!extract_results(step_name, xyz_file)) return false;
        }

        // 构建MECP.x命令
        std::stringstream cmd_ss;
        cmd_ss << "./MECP.x " << step_name;
        cmd_ss << " --tde " << std::scientific << m_input.control.tde;
        cmd_ss << " --tgmax " << std::scientific << m_input.control.tgmax; 
        cmd_ss << " --tgrms " << std::scientific << m_input.control.tgrms;
        cmd_ss << " --tdxmax " << std::scientific << m_input.control.tdxmax;
        cmd_ss << " --tdxrms " << std::scientific << m_input.control.tdxrms;
        cmd_ss << " --stpmx " << std::scientific << m_input.control.stpmx;

        std::string mecp_cmd = cmd_ss.str();

        if (m_debug) {
            std::cout << "Debug: MECP command: " << mecp_cmd << std::endl;
        }

        return execute_command(mecp_cmd, m_debug);
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
            replace_all(tmpl, "[nstep]", std::to_string(step));
            replace_all(tmpl, "[tmpdir]", m_tmpdir);
            if (step == 1) {
                tmpl = std::regex_replace(tmpl, std::regex(R"(\[.*?\])"), "");
            } else {
                tmpl = std::regex_replace(tmpl, std::regex(R"(\[(.*?)\])"), "$1");
            }
            return tmpl;
        };
        std::string input1 = process_template(m_input.inp_tmplt1);
        
        if (m_is_external) {
            if (!write_file(step_name + "_state1.sh", input1)) return false;
            // Only create state2 input if InpTmplt2 is provided
            if (has_inp_tmplt2()) {
                std::string input2 = process_template(m_input.inp_tmplt2);
                if (!write_file(step_name + "_state2.sh", input2)) return false;
            }
        } else {
            std::string suffix = m_config.config["main"]["INPUT_SUFFIX"];
            if (!write_file(step_name + "_state1." + suffix, input1)) return false;
            // Only create state2 input if InpTmplt2 is provided
            if (has_inp_tmplt2()) {
                std::string input2 = process_template(m_input.inp_tmplt2);
                if (!write_file(step_name + "_state2." + suffix, input2)) return false;
            }
        }
        return true;
    }
    
    bool run_qm_calculations(const std::string& step_name) {
        if (m_is_external) {
            if (!execute_command("bash " + step_name + "_state1.sh", m_debug)) return false;
            // Only run state2 if InpTmplt2 is provided
            if (has_inp_tmplt2()) {
                if (!execute_command("bash " + step_name + "_state2.sh", m_debug)) return false;
            } else if (m_debug) {
                std::cout << "Debug: Skipping state2 calculation - InpTmplt2 not provided" << std::endl;
            }
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
            // Only run state2 if InpTmplt2 is provided
            if (has_inp_tmplt2()) {
                if (!run_single_state(2)) return false;
            } else if (m_debug) {
                std::cout << "Debug: Skipping state2 calculation - InpTmplt2 not provided" << std::endl;
            }
        }
        return true;
    }
    
    bool extract_results(const std::string& step_name, const std::string& xyz_file) {
         XYZ current_geom;
        if (!current_geom.read(xyz_file)) return false;
        std::string output_suffix;
        if (m_is_external) {
            // In external mode, try to get OUTPUT_SUFFIX from config, otherwise use default
            if (m_config.config.count("main") > 0 && m_config.config["main"].count("OUTPUT_SUFFIX") > 0) {
                output_suffix = m_config.config["main"]["OUTPUT_SUFFIX"];
            } else {
                output_suffix = "gms"; // Default for external mode
                if (m_debug) {
                    std::cout << "Debug: Using default OUTPUT_SUFFIX 'gms' for external mode" << std::endl;
                }
            }
        } else {
            output_suffix = m_config.config["main"]["OUTPUT_SUFFIX"];
        }
auto extract_single_state = [&](int state_num) {
            std::string state_type = m_input.grp_tmplt.at(state_num == 1 ? "state1" : "state2");
            const auto& rules = m_config.config.at(state_type);
            std::string output_file = step_name + "_state" + std::to_string(state_num) + "." + output_suffix;
            std::string content = read_file(output_file);
            if (content.empty()) return false;
            
            // Extract energy with E.LoacteCount support
            std::smatch e_match;
            std::regex e_regex(rules.at("E"));
            
            // Get E.LoacteCount parameter, default to 1 (first occurrence)
            int e_locate_count = 1;
            if (rules.count("E.LoacteCount") > 0) {
                try {
                    e_locate_count = std::stoi(rules.at("E.LoacteCount"));
                } catch (const std::exception& e) {
                    if (m_debug) {
                        std::cerr << "Debug: Invalid E.LoacteCount, using default value 1" << std::endl;
                    }
                    e_locate_count = 1;
                }
            }
            
            // Find the Nth occurrence of energy pattern
            std::string::const_iterator searchStart = content.cbegin();
            bool found_energy = false;
            for (int i = 0; i < e_locate_count; ++i) {
                if (std::regex_search(searchStart, content.cend(), e_match, e_regex)) {
                    if (i == e_locate_count - 1) {
                        found_energy = true;
                        break;
                    }
                    searchStart = e_match.suffix().first;
                } else {
                    break;
                }
            }
            
            if (!found_energy || e_match.size() < 2) {
                std::cerr << "Error: Could not extract energy from " << output_file 
                          << " (occurrence #" << e_locate_count << ")" << std::endl;
                return false;
            }
            
            std::string energy = e_match[1].str();
            
            if (m_debug) {
                std::cout << "Debug: Found energy (occurrence #" << e_locate_count 
                          << "): " << energy << std::endl;
            }
            
            std::string gard_locate = rules.at("GARD.Loacte");
            
            // Get GARD.LoacteCount parameter, default to 1 (first occurrence)
            int locate_count = 1;
            if (rules.count("GARD.LoacteCount") > 0) {
                try {
                    locate_count = std::stoi(rules.at("GARD.LoacteCount"));
                } catch (const std::exception& e) {
                    if (m_debug) {
                        std::cerr << "Debug: Invalid GARD.LoacteCount, using default value 1" << std::endl;
                    }
                    locate_count = 1;
                }
            }
            
            // Find the Nth occurrence of gard_locate
            size_t gard_pos = std::string::npos;
            size_t search_pos = 0;
            for (int i = 0; i < locate_count; ++i) {
                gard_pos = content.find(gard_locate, search_pos);
                if (gard_pos == std::string::npos) {
                    break;
                }
                search_pos = gard_pos + 1; // Move past current match for next search
            }
            
            if (gard_pos == std::string::npos) {
                std::cerr << "Error: Could not locate gradient block in " << output_file << std::endl;
                if (m_debug) {
                    std::cerr << "Debug: Looking for gradient marker (occurrence #" << locate_count << "): '" << gard_locate << "'" << std::endl;
                }
                return false;
            }
            if (m_debug) {
                std::cout << "Debug: Found gradient block marker at position " << gard_pos << std::endl;
                std::cout << "Debug: Gradient marker: '" << gard_locate << "'" << std::endl;
            }
            std::stringstream content_stream(content.substr(gard_pos));
            std::string line;
            std::getline(content_stream, line);
            if (m_debug) {
                std::cout << "Debug: First line after marker: '" << line << "'" << std::endl;
            }
            int n_skip = std::stoi(rules.at("GARD.NLineSkip"));
            if (m_debug) {
                std::cout << "Debug: Skipping " << n_skip << " lines" << std::endl;
            }
            for (int i = 0; i < n_skip; ++i) {
                std::getline(content_stream, line);
                if (m_debug) {
                    std::cout << "Debug: Skipped line " << (i+1) << ": '" << line << "'" << std::endl;
                }
            }
            std::vector<int> cols;
            std::stringstream ss(rules.at("GARD.TargetColumns"));
            std::string col_str;
            while(std::getline(ss, col_str, ',')) cols.push_back(std::stoi(col_str));
            std::string end_by = rules.at("GARD.EndBy");
            
            // Get GRAD.Type parameter to determine if we need to negate gradients
            // Default is "force" (no negation) for backward compatibility
            std::string grad_type = "force";
            if (rules.count("GRAD.Type") > 0) {
                grad_type = rules.at("GRAD.Type");
                std::transform(grad_type.begin(), grad_type.end(), grad_type.begin(),
                              [](unsigned char c){ return std::tolower(c); });
            }
            
            // Determine if gradients need to be negated
            bool negate_gradients = (grad_type == "grad");
            double gradient_factor = negate_gradients ? -1.0 : 1.0;
            
            if (m_debug) {
                std::cout << "Debug: GRAD.Type = '" << grad_type << "', gradient_factor = " 
                          << gradient_factor << std::endl;
            }
            
            std::stringstream grad_file_content;
            grad_file_content << current_geom.atom_count << "\n";
            try {
                grad_file_content << std::fixed << std::setprecision(10) << std::stod(energy) << "\n";
            } catch (const std::exception& e) {
                std::cerr << "Error: Failed to convert energy string to double." << std::endl;
                std::cerr << "Energy string: '" << energy << "'" << std::endl;
                std::cerr << "Exception: " << e.what() << std::endl;
                return false;
            }
            for (int i = 0; i < current_geom.atom_count; ++i) {
                std::getline(content_stream, line);
                if (m_debug) {
                    std::cout << "Debug: Reading gradient line " << (i+1) << " for atom " << current_geom.symbols[i] << ": '" << line << "'" << std::endl;
                }
                if ((!end_by.empty() && line.find(end_by) != std::string::npos) || (end_by.empty() && line.find_first_not_of(" \t\n\r") == std::string::npos)) {
                    std::cerr << "Error: Premature end of gradient block." << std::endl;
                    return false;
                }
                std::stringstream line_ss(line);
                std::vector<std::string> tokens;
                std::string token;
                while(line_ss >> token) tokens.push_back(token);
                grad_file_content << current_geom.symbols[i] << "    ";
                
                try {
                    double grad_x = std::stod(tokens.at(cols[0]-1)) * gradient_factor;
                    double grad_y = std::stod(tokens.at(cols[1]-1)) * gradient_factor;
                    double grad_z = std::stod(tokens.at(cols[2]-1)) * gradient_factor;

                    grad_file_content << std::fixed << std::setprecision(10) 
                                     << grad_x << "    " << grad_y << "    " << grad_z << "\n";
                } catch (const std::exception& e) {
                    std::cerr << "Error: Failed to parse gradient for atom " << (i+1) 
                              << " in " << output_file << std::endl;
                    if (m_debug) {
                        std::cerr << "Debug: " << e.what() << " - Line: '" << line << "'" << std::endl;
                    }
                    return false;
                }
            }
            return write_file(step_name + ".grad" + std::to_string(state_num), grad_file_content.str());
        };
        if (!extract_single_state(1)) return false;
        // Always extract second state results, regardless of whether InpTmplt2 was provided
        if (!extract_single_state(2)) return false;
        return true;
    }

    bool check_convergence() {
        if (!fs::exists("convg.tmp")) {
             std::cerr << "Warning: convg.tmp not found. Assuming not converged." << std::endl;
            return false;
        }

        // Read only the first line to get the status
        std::ifstream file("convg.tmp");
        if (!file.is_open()) {
            std::cerr << "Warning: Cannot open convg.tmp. Assuming not converged." << std::endl;
            return false;
        }

        std::string status;
        std::getline(file, status);
        file.close();

        // Trim whitespace from status string
        status.erase(0, status.find_first_not_of(" \t\n\r"));
        status.erase(status.find_last_not_of(" \t\n\r") + 1);

        if (m_debug) {
            std::cout << "Debug: Convergence status from first line: '" << status << "'" << std::endl;
        }

        if (status == "CONVERGED") {
            if (!fs::exists("final.xyz")) {
                std::cerr << "Error: Convergence indicated, but final.xyz was not created." << std::endl;
                return false; // Treat as error/not converged
            }
            append_to_trajectory("final.xyz");
            return true;
        } else {
            // Only append new.xyz if it exists
            if (fs::exists("new.xyz")) {
                append_to_trajectory("new.xyz");
            } else if (m_debug) {
                std::cout << "Debug: new.xyz not found for trajectory append" << std::endl;
            }
            return false;
        }
    }
    
    std::string parse_convergence_info() {
        if (!fs::exists("convg.tmp")) {
            return "Convergence info not available";
        }
        
        std::string content = read_file("convg.tmp");
        if (content.empty()) {
            return "Convergence info not available";
        }
        
        // Parse the convergence information
        std::map<std::string, std::string> conv_data;
        std::stringstream ss(content);
        std::string line;
        
        while (std::getline(ss, line)) {
            // Remove leading/trailing whitespace
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);
            
            if (line.empty()) continue;
            
            // Parse different types of lines
            if (line.find("Energy_Gap:") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string value = line.substr(pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    conv_data["Energy_Gap"] = value;
                }
            } else if (line.find("Max_Gradient:") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string value = line.substr(pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    // Extract just the current value (before any parentheses)
                    size_t paren_pos = value.find('(');
                    if (paren_pos != std::string::npos) {
                        value = value.substr(0, paren_pos);
                        value.erase(value.find_last_not_of(" \t") + 1);
                    }
                    conv_data["Max_Gradient"] = value;
                }
            } else if (line.find("RMS_Gradient:") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string value = line.substr(pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    size_t paren_pos = value.find('(');
                    if (paren_pos != std::string::npos) {
                        value = value.substr(0, paren_pos);
                        value.erase(value.find_last_not_of(" \t") + 1);
                    }
                    conv_data["RMS_Gradient"] = value;
                }
            } else if (line.find("Max_Displacement:") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string value = line.substr(pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    size_t paren_pos = value.find('(');
                    if (paren_pos != std::string::npos) {
                        value = value.substr(0, paren_pos);
                        value.erase(value.find_last_not_of(" \t") + 1);
                    }
                    conv_data["Max_Displacement"] = value;
                }
            } else if (line.find("RMS_Displacement:") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string value = line.substr(pos + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    size_t paren_pos = value.find('(');
                    if (paren_pos != std::string::npos) {
                        value = value.substr(0, paren_pos);
                        value.erase(value.find_last_not_of(" \t") + 1);
                    }
                    conv_data["RMS_Displacement"] = value;
                }
            }
        }
        
        // Format the convergence information
        std::stringstream formatted;
        formatted << "E=" << (conv_data.count("Energy_Gap") ? conv_data["Energy_Gap"] : "N/A");
        formatted << " MaxF=" << (conv_data.count("Max_Gradient") ? conv_data["Max_Gradient"] : "N/A");
        formatted << " RMSF=" << (conv_data.count("RMS_Gradient") ? conv_data["RMS_Gradient"] : "N/A");
        formatted << " MaxD=" << (conv_data.count("Max_Displacement") ? conv_data["Max_Displacement"] : "N/A");
        formatted << " RMSD=" << (conv_data.count("RMS_Displacement") ? conv_data["RMS_Displacement"] : "N/A");
        
        return formatted.str();
    }

    void append_to_trajectory(const fs::path& xyz_path) {
        if (!fs::exists(xyz_path)) {
            std::cerr << "Warning: Cannot append to trajectory, file not found: " << xyz_path << std::endl;
            return;
        }
        
        // Read the xyz file
        std::ifstream xyz_file(xyz_path);
        if (!xyz_file.is_open()) {
            std::cerr << "Warning: Cannot open xyz file for trajectory: " << xyz_path << std::endl;
            return;
        }
        
        // Read and modify the xyz content
        std::stringstream modified_content;
        std::string line;
        int line_count = 0;
        
        while (std::getline(xyz_file, line)) {
            line_count++;
            if (line_count == 2) {  // Second line is the comment line in xyz format
                // Replace comment with convergence information
                std::string conv_info = parse_convergence_info();
                modified_content << conv_info << "\n";
                if (m_debug) {
                    std::cout << "Debug: Replaced comment line with: " << conv_info << std::endl;
                }
            } else {
                modified_content << line << "\n";
            }
        }
        xyz_file.close();
        
        // Append to trajectory file
        std::ofstream traj_file("opt.trj", std::ios::app);
        if (traj_file.is_open()) {
            traj_file << modified_content.str();
            if (m_debug) {
                std::cout << "Debug: Appended modified structure to trajectory file" << std::endl;
            }
        } else {
            std::cerr << "Warning: Cannot open trajectory file for writing" << std::endl;
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
        std::cerr << "    restart=true" << std::endl;
        std::cerr << "  end" << std::endl;
        std::cerr << "  %InpTmplt1 ... end" << std::endl;
        std::cerr << "  %InpTmplt2 ... end (optional)" << std::endl;
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