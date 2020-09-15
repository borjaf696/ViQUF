//
// Created by borja on 31/3/20.
//

#ifndef VIADBG_TRIAL_EXTRA_H
#define VIADBG_TRIAL_EXTRA_H

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
    string t_data = "virus";
    double missmatches;
    size_t genome_size = 0;
    double num_unique_kmers = 0;
    size_t accumulative_h = 0;
    bool full_info = false, metagenomic = false,
            postProcess = false, remove_duplicates = true, polish = true, gfa = false, debug = false, greedy = false;
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
};
/*
 * Maths operations
 */
struct Maths{

    template<typename T>
    static T my_round(T value, bool up = false)
    {
        float result = value, first_val = 1, jump = 10;
        while (result > MARGIN)
        {
            first_val *= jump;
            result = value / first_val;
        }
        return ((!up)?floor(value / first_val):ceil(value/first_val)) * first_val;
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
        return container[ceil((float) container.size() * 0.5)-1];
    }
};


struct Common {
    static void display_unitig(const vector<string> & sequence_map, size_t place)
    {
        bool reverse = place >= sequence_map.size();
        string seq = sequence_map.at(reverse?place - sequence_map.size():place);
        seq = ((reverse)?Sequence(seq.c_str()).getRevcomp():seq);
        cout<<" "<<seq<<" ";
    }

    static string return_unitig(const vector<string> & sequence_map, size_t place)
    {
        bool reverse = (place >= sequence_map.size());
        string seq = sequence_map.at(reverse?place - sequence_map.size():place);
        return ((reverse)?Sequence(seq.c_str()).getRevcomp():seq);
    }

    static size_t return_index(const vector<string> & sequence_map, size_t place)
    {
        bool reverse = (place >= sequence_map.size());
        return (place >= sequence_map.size())?place - sequence_map.size():place;
    }
};

struct System {
};
#endif //VIADBG_TRIAL_EXTRA_H
