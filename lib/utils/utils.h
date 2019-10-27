#ifndef   UTILS_H
#define   UTILS_H
    #include <vector>
    #include <set>
    #include <string>
    #include <fstream>
    #include <sstream>
    #include <stdio.h>
    #include <math.h>
    #include <fstream>
    #include <iostream>
    #include <sys/stat.h>
    #include <boost/config.hpp>
    #include <boost/filesystem.hpp>
    #include <boost/program_options.hpp>
    #include <boost/program_options/detail/config_file.hpp>
    #include <boost/program_options/parsers.hpp>
    #include "sdsl/construct.hpp"
    #include "sdsl/suffix_arrays.hpp"
    #include "sdsl/bit_vectors.hpp"
    #include "sdsl/csa_bitcompressed.hpp"
    #include "sdsl/lcp_bitcompressed.hpp"
    #include "OptionPrinter.hpp"

    using namespace std;

    template <class t_array>
    void print_array(const t_array &array, string name = ""){
        for(uint64_t i = 0; i < array.size(); ++i){
            cout << name << "[" << i << "]=" << array[i] << endl;
        }
    }

    template <class t_array>
    void check_lcp (const t_array &array, size_t threshold)
    {
        for (size_t i = 0; i < array.size(); ++i)
        {
            if (array[i] < threshold)
                cout << "[" << i << "]=" <<array[i] <<endl;
        }
    }
    /*
     * Progress
     */
    struct Progress
    {
        static void update(size_t num_actual){
            if (Progress::get().show) {
                int val=(int) (((float) num_actual) / ((float) Progress::get().size_total) * 100);
                if (val % 10 == 0)
                    if (!Progress::get().f[val/10]) {
                        Progress::get().f[val/10] = true;
                        show_progress(num_actual);
                    }
            }
        }

        static Progress& get(){
            static Progress progress;
            return progress;
        }
        std::vector<bool> f = std::vector<bool>(10,false);
        size_t size_total;
        bool show = false;
    private:
        static void show_progress(size_t num_actual){
            size_t val = (size_t)((float)num_actual/
                                  (float)Progress::get().size_total * 100);
            std::cout << val << ((val==100)?"%\n":"% ") << std::flush;
            if (val == 100) {
                _prepare_next();
                /*
                 * Stop showing results
                 */
                Progress::get().show = false;
            }
        }
        static void _prepare_next(){
            Progress::get().f = std::vector<bool>(10,false);
        }
    };
    /*
     * System operations
     */
    struct System{
        /*
         * System utils
         */
        static void removeFile(string filename)
        {
            remove(filename.c_str());
        }
        static void removeDir(string dirname)
        {
            namespace fs = boost::filesystem;
            fs::remove_all(dirname.c_str());
        }
        static void createdir(string dirname)
        {
            int check = mkdir(dirname.c_str(),0777);
            if (check!=-1)
                cout << "Directory created "<<dirname<<endl;
            else {
                cout << "Unable to create directory "<<dirname<<endl;
            }
        }
        static void changeName(string originalName, string newName)
        {
            rename(originalName.c_str(), newName.c_str());
        }
        static string findExtension(string s)
        {
            return s.substr(s.find_last_of(".") + 1);
        }
        static bool exist(string && s)
        {
            struct stat buffer;
            return (stat(s.c_str(), &buffer) == 0);
        }
        static set<string> getAllFaFqFiles(string path, string chain = "", bool recursive = false)
        {
            namespace fs = boost::filesystem;
            set<string> output;
            fs::path fs_path(path);
            if (fs::is_regular_file(fs_path))
                output.emplace(path);
            else
            {
                for (auto & p : fs::directory_iterator(path))
                {
                    if (fs::is_regular_file(p))
                    {
                        bool add = true;
                        if (chain != "" && (chain != p.path().string()))
                            add = false;
                        if (add)
                        {
                            ostringstream oss;
                            oss << p;
                            string converted_path = oss.str().substr(1, oss.str().size() - 2);
                            string extension = converted_path.substr(converted_path.rfind('.'));
                            if (extension == ".fastq" || extension == ".fq" || extension == ".fasta" || extension == ".fa")
                                output.emplace(converted_path);
                        }
                    }else if (recursive)
                    {
                        set<string> files = getAllFaFqFiles(p.path().string(),chain,recursive);
                        output.insert(files.begin(),files.end());
                    }
                }
            }
            return output;
        }
        static void execute(std::string instruction)
        {
            if (system(instruction.c_str()))
            {
                cout << "Fail on: "<<instruction<<"\n";
                exit(1);
            }
        }
        static std::string appendFiles(std::vector<string> files, std::string newFile)
        {
            std::string instruction = "cat ";
            for (auto s: files)
                instruction += s+" ";
            instruction += ">"+newFile;
            cout << "Instruction: "<<instruction<<endl;
            if (system(instruction.c_str())){
                cout << "Problem executing: "<<instruction<<"\n";
                cout << "Exiting\n";
                exit(1);
            }
            return newFile;
        }
        template<typename T>
        static void write_histogram(std::unordered_map<T,std::vector<size_t>> histogram, std::string file_name)
        {
            ofstream file;
            file.open(file_name);
            for (auto k: histogram)
            {
                file << (k.first.str()+"\n");
                for (auto t: k.second)
                    file << (std::to_string(t)+" ");
                file << endl;
            }
            file.close();
        }
    };
    /*
     * Basic Operations
     */
    struct Basics
    {
        static float mean(vector<int> contigsLengthVector)
        {
            int count = 0;
            for (auto v:contigsLengthVector)
                count += v;
            return count / contigsLengthVector.size();
        }
        static float standardDeviation(vector<int> contigsLengthVector, float mean)
        {
            float var = 0;
            for (auto v: contigsLengthVector)
                var += (v-mean)*(v-mean);
            var /= contigsLengthVector.size();
            return sqrt(var);
        }
    };
    /*
     * Bioinformatics
     */
    struct Bio
    {
        static bool pairedEnd(string A, string B)
        {
            string delimiter = ".";
            string token, token2;
            size_t posA, posB;
            while (((posA = A.find(delimiter)) != string::npos) && ((posB = B.find(delimiter)) != string::npos)) {
                token = A.substr(0, posA);
                token2 = B.substr(0, posB);
                cout << token << " "<<token2<<endl;
                if (token == token2) {
                    A.erase(0, posA + delimiter.length());
                    B.erase(0, posB + delimiter.length());
                }else{
                    A.replace(A.size(),1,1,'2');
                    cout << A <<" "<<B<<endl;
                    if (A == B)
                        return true;
                    return false;
                }
            }
            return false;
        }
        static bool isFasta(string s)
        {
            return (System::findExtension(s) == "fa" || System::findExtension(s) == "fasta");
        }
        template<typename T>
        static void dskCounter(unordered_map<T, size_t> & T_maps, size_t & id, string && file)
        {
            size_t total_lines = 0, num_line = 0;
            ifstream infile(file);
            for (string line;getline(infile,line);)
                total_lines++;
            infile.close();
            /*
             * Progress
             */
            {
                Progress::get().show = true;
                Progress::get().size_total = total_lines;
            }
            infile.open(file);
            for( string line; getline( infile, line ); )
            {
                if (line.back() == '\n' or line.back() == '\r')
                {
                    line.pop_back();
                }
                string delimiter = " ", token;
                size_t pos, times = 0;
                T y;
                while ((pos = line.find(delimiter)) != string::npos)
                {
                    token = line.substr(0,pos);
                    if (times++ == 0)
                    {
                        y = T(token);
                        break;
                    }
                    line.erase(0, pos + delimiter.length());
                }
                if (T_maps.find(y) == T_maps.end())
                    T_maps[y] = id++;
                Progress::update(num_line++);
            }
            Progress::update(total_lines);
            infile.close();
        }
        template<typename T, typename G>
        static size_t dskParser(unordered_map<G, unordered_map<T,size_t>> & map,
                unordered_map<T, size_t> & T_maps, string && file, G environment)
        {
            size_t num_line = 0, total_lines = 0, max_freq = 0;
            ifstream infile(file);
            for (string line;getline(infile,line);)
                total_lines++;
            infile.close();
            /*
             * Progress
             */
            {
                Progress::get().show = true;
                Progress::get().size_total = total_lines;
            }
            infile.open(file);
            for( string line; getline( infile, line ); )
            {
                if (line.back() == '\n' or line.back() == '\r')
                {
                    line.pop_back();
                }
                string delimiter = " ", token;
                size_t pos, times = 0;
                T y;
                while ((pos = line.find(delimiter)) != string::npos)
                {
                    token = line.substr(0,pos);
                    if (times++ == 0)
                        y = T(token);
                    line.erase(0, pos + delimiter.length());
                }
                times = atoi(line.c_str());
                max_freq = (max_freq < times)?times:max_freq;
                map[environment][y] = times;
                Progress::update(num_line++);
            }
            Progress::update(total_lines);
            infile.close();
            return max_freq;
        }
        template<typename T>
        static void fastaWriter(vector<T> seqVector, string outputFile)
        {
            outputFile+=".fasta";
            ofstream singleFile(outputFile, ios::app);
            size_t num_kmer = 0;
            for (auto s:seqVector)
            {
                singleFile << '>'+to_string(num_kmer++)<<endl;
                singleFile << s.str()<<endl;
            }
        }

    };
#endif