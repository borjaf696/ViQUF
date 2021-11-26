//
// Created by borja on 31/3/20.
//

#ifndef VIADBG_TRIAL_EXTRA_H
#define VIADBG_TRIAL_EXTRA_H
#include <filesystem>
#include <sys/stat.h>
// For viral data it must be higher (100) for meta try lower
#define MARGIN 100
/*
 * Parameters
 */
struct Parameters
{
    static Parameters& get(){
        static Parameters param;
        return param;
    }
    void set(size_t param, std::string param_read)
    {
        if (param_read == "KmerSize")
            Parameters::get().kmerSize = param;
        if (param_read == "Accumulative")
            Parameters::get().accumulative_h = param;
        if (param_read == "missmatches")
            Parameters::get().missmatches = param;
        if (param_read == "numKmers")
            Parameters::get().num_unique_kmers = param;
        if (param_read == "genome_size")
            Parameters::get().genome_size = param;
        if (param_read == "full_info")
            Parameters::get().full_info = param;
        if (param_read == "metagenomic")
            Parameters::get().metagenomic = param;
        if (param_read == "remove_duplicates")
            Parameters::get().remove_duplicates = param;
        if (param_read == "numThreads")
            Parameters::get().numThreads = param;
        if (param_read == "polish")
            Parameters::get().polish = param;
        if (param_read == "postProcess")
            Parameters::get().postProcess = param;
        if (param_read == "gfa")
            Parameters::get().gfa = param;
        if (param_read == "debug")
            Parameters::get().debug = param;
        if (param_read == "greedy")
            Parameters::get().greedy = param;
        if (param_read == "t_data")
            Parameters::get().t_data = param;
    }

    void check_cmd_line(int argc, char* argv[])
    {
        for (size_t i = 0; i < argc; i++)
        {
            if (strcmp(argv[i],"--debug") == 0)
                Parameters::get().debug = true;
            if (strcmp(argv[i],"--greedy") == 0)
                Parameters::get().greedy = true;
            if (strcmp(argv[i],"--virus") == 0)
                Parameters::get().t_data = "virus";
            if (strcmp(argv[i],"--Illumina") == 0)
                Parameters::get().t_data = "Illumina";
            if (strcmp(argv[i],"--Amplicons") == 0)
                Parameters::get().t_data = "Amplicons";
            if (strcmp(argv[i],"--TGS") == 0)
                Parameters::get().t_data = "TGS";
            if (strcmp(argv[i],"--amplicon") == 0)
                Parameters::get().partial_assembly = true;
        }
        Parameters::get().kmerSize = atoi(argv[2]);
    }

    void print_info()
    {
        cout << "################ Execution Information #################"<<endl;
        cout << "Data type: "<<t_data<<endl;
        cout << "Kmer size: "<<kmerSize<<endl;
        cout << "Metagenomic parameters: "<<metagenomic<<endl;
        cout << "Debugging: "<<debug<<endl;
        cout << "Greedy extension: "<<greedy<<endl;
        cout << "########################################################"<<endl;
    }
    string t_data = "virus", paired_end_reads = "";
    double missmatches;
    size_t genome_size = 0;
    double num_unique_kmers = 0;
    size_t accumulative_h = 0;
    bool full_info = false, metagenomic = false,
            postProcess = false, remove_duplicates = true, 
            polish = true, gfa = false, debug = false, greedy = false, partial_assembly = false;
    size_t kmerSize;
    size_t numThreads;
    bool show = false;
};
/*
 * Containers operations
 */
struct Op_Container{
    template<typename T>
    static void join_containers(vector<T> & src, const vector<T> & target)
    {
        for (auto s: target)
            src.push_back(s);
    }
    template<typename T>
    static unordered_set<T> to_set(const vector<T> & src )
    {
        unordered_set<T> out;
        for (auto s: src)
            out.emplace(s);
        return out;
    }
    template<typename T>
    static bool in(const vector<T> & v1, T el)
    {
        return (std::find(v1.begin(), v1.end(), el) != v1.end());
    }
    template<typename T>
    static pair<size_t, float> percentage_of_similarity(const vector<T> & v1, const vector<T> & v2)
    {
        float matches = 0;
        for (auto v: v1)
            matches += (Op_Container::in(v2, v))?1:0;
        // 1 if v2 in v1 0 otherwise
        return (v1.size() > v2.size())?std::pair<size_t,float>(1, matches/v2.size()):std::pair<size_t, float>(0, matches/v1.size());
    }
};
/*
 * Maths operations
 */
struct Maths{
    static float function_val( unordered_map<size_t,size_t> distribution, 
        unordered_map<size_t,size_t> survavility,unordered_map<size_t,size_t> length_reads, size_t place, size_t genome_size)
    {
        float term1 = place * (survavility[place] - survavility[genome_size - place]);
        float term2 = 0;
        for (int i = place; i >= 0; --i)
            term2 += (length_reads[i]*i);
        float term3 = (genome_size - place) * (survavility[genome_size - place]);
        return (term1 + term2 + term3);
    }

    static float calculate_readjustment_ratio(unordered_map<size_t,size_t> distribution, 
        unordered_map<size_t,size_t> survavility, unordered_map<size_t,size_t> length_reads, size_t place, size_t genome_size, float average_read_length)
        {
            // Avoid divide by 0
            place = (place == 0)?1:place;
            size_t place_tmp = place;
            if (place > genome_size / 2)
                place_tmp = genome_size - place;
            float val_max = function_val(distribution, survavility,length_reads, (genome_size - average_read_length), genome_size), 
                val = function_val(distribution, survavility,length_reads, place_tmp, genome_size);
            return val_max/val;
        }
    template<typename T>
    static T my_round(T value, bool up = false)
    {
        float result = value, first_val = 10, jump = 10;
        while (result > MARGIN)
        {
            first_val *= jump;
            result = value / first_val;
        }
        return ((!up)?floor(value / first_val):round(value/first_val)) * first_val;
    }

    /*
     * Input must be shorted before trying
     */
    template<typename T>
    static T i_quantile(vector<T> container, float quantile = 0.0, bool next_value = true)
    {
        size_t index = ceil((float) container.size() * quantile);
        T result = container[index];
        if (quantile == 0.0)
            return result;
        while (container[index] == result)
        {
            (next_value)?index++:index--;
        }
        return container[index];
    }
    /*
     * Median
     */
    template<typename T>
    static T median(vector<T> container)
    {
        size_t place = ceil((float) container.size() * 0.5)-1;
        return (container[place] == 0)?(container[place + 1]):container[place];
    }
    template<typename T>
    static T max(vector<T> container)
    {
        T maximum = 0;
        for (auto t:container)
            maximum = (t > maximum)?t:maximum;
        return maximum;
    }
    template<typename T>
    static T min(vector<T> container)
    {
        T minimum = 9999999999;
        for (auto t:container)
            minimum = (t < minimum)?t:minimum;
        return minimum;
    }
};


struct Common {
    static void display_unitig(const vector<string> & sequence_map, size_t place, bool full_unitig_map = false)
    {
        if (!full_unitig_map){
            bool reverse = place >= sequence_map.size();
            string seq = sequence_map.at(reverse?place - sequence_map.size():place);
            char * seq_str = &seq.c_str()[0];
            seq = ((reverse)?Sequence(seq_str).getRevcomp():seq);
            cout<<" "<<seq<<" ";
        } else {
            cout << sequence_map.at(place) << " ";
        }
    }

    static string return_unitig(const vector<string> & sequence_map, size_t place, bool full_unitig_map = false)
    {
        if (!full_unitig_map){
            bool reverse = (place >= sequence_map.size());
            string seq = sequence_map.at(reverse?place - sequence_map.size():place);
            char * seq_str = &(seq.c_str())[0];
            return ((reverse)?Sequence(seq_str).getRevcomp():seq);
        } else {
            return sequence_map.at(place);;
        }
    }

    static size_t return_index(const vector<string> & sequence_map, size_t place, bool full_unitig_map = false)
    {
        if (!full_unitig_map){
            bool reverse = (place >= sequence_map.size());
            return (place >= sequence_map.size())?place - sequence_map.size():place;
        } else {
            return place;
        }
    }
    static size_t rev_comp_index(const size_t total_unitigs, size_t place, bool full_unitig_map = false)
    {
        if (!full_unitig_map)
            return (place >= total_unitigs/2)?place - total_unitigs/2:place + total_unitigs/2;
        return place;
    }

    template<typename T>
    static float sum_vector(const vector<T> & v1)
    {
        float total = 0;
        for (auto v:v1)
            total += (float) v;
        return (float) total;
    }
};
#include <sys/types.h>
#include <dirent.h>
struct OwnSystem {
    static std::vector<string> get_directories(std::string path, std::string extension = "fastq")
    {
        DIR* dirp = opendir(path.c_str());
        struct dirent * dp;
        std::vector<string> files;
        while ((dp = readdir(dirp)) != NULL) {
            if (dp->d_name[0] == '.')
                continue;
            char * p_1 = strtok(dp->d_name,".");
            char * p_2 = strtok(NULL,".");
            std::cout << p_1 << " "<<p_2<<std::endl;
            if (p_2 != NULL){
                if (strcmp(p_2,extension.c_str())==0){
                    std::cout << "Insert: "<<(path+dp->d_name+".fastq")<<std::endl;
                    files.push_back(path+dp->d_name+".fastq");
                }
                std::string extension_fasta = "fasta";
                if (strcmp(p_2,extension_fasta.c_str())==0){
                    std::cout << "Insert: "<<(path+dp->d_name+".fasta")<<std::endl;
                    files.push_back(path+dp->d_name+".fasta");
                }
            }
        }
        closedir(dirp);
        std::sort(files.begin(), files.end());
        return files;
    }
};
#endif //VIADBG_TRIAL_EXTRA_H
