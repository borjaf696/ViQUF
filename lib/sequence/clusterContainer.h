//
// Created by borja on 7/06/19.
//

#ifndef ESPRELA_CLUSTERCONTAINER_H
#define ESPRELA_CLUSTERCONTAINER_H
#include <chrono>
#ifndef UTILS_H
    #include "../utils/utils.h"
#endif
#ifndef  ESPRELA_PROBTABLE_H
    #include "../fs/probTable.h"
#endif
#ifndef SEQUENCE_CONTAINER_H
    #include "../sequence/sequence_container.h"
#endif
#ifndef KMER_H
    #include "../sequence/kmer.h"
#endif
#define MIN_ABUNDANCE 1
#define SAMPLING 1
#define NUM_SAMPLES 15
#define SAMPLE_RATE 0.4

using namespace sdsl;
using namespace std;

typedef size_t envi;

class ClusterContainer{
public:
    ClusterContainer(vector<string> reads, vector<string> f_environments, string directory, string binner, bool clustering, rank_support_v<1> rs_ds):
        _reads(reads),_f_environments(f_environments),_binner(binner),_directory(directory),_clustering(clustering),_rs_ds(rs_ds)
    {
        _doClustering();
        _loadClusters();
    }
    vector<vector<int>> getClusters();
    void classifyClusters();
    void selectKmers()
    {
        _processHybridClusters();
    }
private:
    void _doClustering()
    {
        if (_clustering) {
            cout << "Clustering with "<<_binner<<endl;
            if (_binner == "maxbin")
                _binnerExecution(_directory+"output");
            else if (_binner == "metabat")
                _binnerExecution(_directory+"output");
            else
            {
                cout << "Binner not accepted"<<endl;
                exit(1);
            }
        }else
            cout << "Skipping clustering"<<endl;
    }

    void _loadClusters()
    {
        set<string> cluster_files = System::getAllFaFqFiles(_directory);
        for (auto f: cluster_files)
        {
            cout << "File: "<<f<<endl;
            vector<int> current_cluster;
            ifstream  contigsFile (f, ios::out | ios::app | ios::binary);
            if (contigsFile.is_open())
            {
                string line;
                while ( getline (contigsFile,line) )
                {
                    if ( line[0] == '>'){
                        current_cluster.push_back(atoi(line.substr(1, line.size()).c_str()));
                    }
                }
            }
            _clusters.push_back(current_cluster);
        }
    }

    void _binnerExecution(string output)
    {
        string contig_file = "contigs.fasta", instruction;
        if (_binner == "maxbin")
        {
            instruction = "bash -c \"perl third-party/MaxBin/run_MaxBin.pl -contig " + contig_file + " ";
            instruction += "-reads " + _reads[0] + " ";
            for (size_t i = 1; i < _reads.size(); ++i)
                instruction += "-reads" + to_string(i + 1) + " " + _reads[i] + " ";
            instruction += "-thread 32 -out " + output + "\"";
            cout << "Instruction: " << instruction << endl;
        }
        if (_binner == "metabat")
        {
            string contig_file = "contigs.fasta";
            string instruction = "export PATH=$PATH:~/metabat/bin/;bash -c \"metabat2 -i "+contig_file+" -o "+output+"\"";
            cout << "Instruction: "<<instruction<<endl;
        }
        System::execute(instruction);
    }

    void _processHybridClusters()
    {
        set<string> cluster_files = System::getAllFaFqFiles(_directory);
        vector<string> cluster_files_vect;
        cluster_files_vect.insert(cluster_files_vect.end(), cluster_files.begin(), cluster_files.end());
        for (size_t i = 0; i < _hybridClusters.size(); i++)
            if (_hybridClusters[i])
                _selectSignificantKmersForClusters(cluster_files_vect[i]);
    }

    void _selectSignificantKmersForClusters(string clusterFile)
    {
        cout << "Process: "<<clusterFile<<endl;
        //FSContext fsContext(clusterFile);
    }

    vector<string> _reads,_f_environments;
    const string _binner, _directory;
    bool _clustering;
    rank_support_v<1> _rs_ds;
    vector<vector<int>> _clusters;
    vector<bool> _hybridClusters;
};

class SimpleEnvironment{
public:
    ~SimpleEnvironment()
    {
        delete _fsContext;
        delete _probTable;
    }
    SimpleEnvironment(vector<string> reads, vector<string> environments):_reads(reads),
                                            _environments(environments), _reads_processed(0), _last_kmer(0)
    {
        envi num_env = 0;
        string curr_env = environments[0];
        _samples_map[curr_env] = num_env;
        for(auto e:_environments)
        {
            if (e != curr_env) {
                curr_env = e;
                _samples_map[curr_env] = ++num_env;
            }
        }
        /*
         * Correccion lecturas
         */
        vector<string> cumulative_reads, corrected_files;
        string newFile;
        num_env = _samples_map[environments[0]];
        bool fasta = Bio::isFasta(_reads[0]);
        /*
         * Timing
         */
        auto start_corr = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < _environments.size(); ++i)
        {
            if (num_env == _samples_map[_environments[i]])
                cumulative_reads.push_back(_reads[i]);
            else{
                for (auto r: cumulative_reads)
                    cout << "File: "<<r<<endl;
                newFile = System::appendFiles(cumulative_reads, "output_file.fasta");
                corrected_files.push_back(_correctReads(newFile, "corrected_"+to_string(num_env)+((fasta)?".fa":".fq")));
                num_env = _samples_map[_environments[i]];
                cumulative_reads.clear();
                cumulative_reads.push_back(_reads[i]);
            }
        }
        newFile = System::appendFiles(cumulative_reads, "output_file.fasta");
        corrected_files.push_back(_correctReads(newFile, "corrected_"+to_string(num_env)+((fasta)?".fa":".fq")));
        auto finish_corr = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_corr = finish_corr - start_corr;
        cout << "Elapsed Correction: " << elapsed_corr.count() << " s\n";
        /*
         * Cuenta de k-mers
         */
        newFile = System::appendFiles(corrected_files, "all_reads.fasta");
        _launchCounter("all_reads.fasta", 0, true);
        _vector_kmer_tmp = vector<Kmer>(_kmers_map.size());
        for (auto k:_kmers_map) {
            _vector_kmer_tmp[k.second] = k.first;
        }
        cout << "Stats: "<<_kmers_map.size()<<" "<<NUM_SAMPLES*(num_env+1)<<endl;
        /*
         * Feature selection context
         */
        _fsContext = new FSContext(_kmers_map.size(), NUM_SAMPLES*(num_env+1));
        _probTable = new ProbTable(_kmers_map.size());
        /*
         * Sampling
         */
        num_env = _samples_map[_environments[0]];
        size_t num_sample = 0;
        /*
         * Ojo con muestras del mismo environment pero distintas samples
         */
        auto start_samp = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < corrected_files.size(); ++i)
        {
            cout << "NumEnv: "<<num_env<<endl;
            _initialization(corrected_files[i], num_sample, num_env);
            num_env = i+1;
        }
        _probTable->calculateProbTable(*(_fsContext));
        //_fsContext->showInfo();
        cout << "End Basics"<<endl;
        auto finish_samp = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_samp = finish_samp - start_samp;
        cout << "Elapsed Sampling: " << elapsed_samp.count() << " s\n";
        /*
         * Fast-mRMR
         */
        auto start_fast = std::chrono::high_resolution_clock::now();
        cout << "Initializing fast-mRMR"<<endl;
        fastmRMR fm(_probTable,_fsContext, num_env+1);
        cout << "End initialization fast-mRMR"<<endl;
        vector<size_t> fselection = fm.getFeatures(0.8);
        vector<Kmer> kmers_selected;
        cout << "Features 1: "<<endl;
        for (auto f:fselection) {
            cout << f << ", ";
            kmers_selected.push_back(_vector_kmer_tmp[f]);
        }
        _reportKmers(kmers_selected);
        auto finish_fast = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_fast = finish_fast - start_fast;
        cout << "Elapsed fast-mRMR: " << elapsed_fast.count() << " s\n";
        cout << endl;
    }
private:
    /*
     * Report kmers
     */
    void _reportKmers(vector<Kmer> kmersSelected)
    {
        string outputFile = "selectedKmers";
        System::removeFile(outputFile+".fasta");
        Bio::fastaWriter(kmersSelected, outputFile);
    }
    /*
     * First approach:
     *  Contar por sample
     *  Crear matriz
     *  Normalizar y discretizar
     */
    void _initialization(string newFile, size_t & num_sample, size_t num_env)
    {
        SequenceContainer sc_temp;
        sc_temp.loadFromFile(newFile, false, num_env);
        size_t max_freq_kmer = 0;
        if (SAMPLING)
        {
            /*
             * Sampling sequence container
             */
            System::removeDir("sample/");
            cout << "Starting sampling"<<endl;
            string path = sc_temp.sampling(SAMPLE_RATE, NUM_SAMPLES, to_string(num_env));
            cout << "End sampling" << endl;
            set <string> sub_samples = System::getAllFaFqFiles(path);
            for (auto sample: sub_samples) {
                cout << "Counting sample: "<<sample<<endl;
                max_freq_kmer = _launchCounter(sample, num_sample);
                cout << "Matrix building starts for sample: "<<num_sample<<endl;
                _samples_kmers_map[num_sample] = max_freq_kmer;
                _checkSamples();
                cout << "Matrix building ends for sample: "<<num_sample++<< endl;
                _envi_kmers_map.clear();
            }
        }else {
            max_freq_kmer = _launchCounter(newFile, num_env);
            _samples_kmers_map[num_env] = max_freq_kmer;
        }
    }
    string _correctReads(string file, string output)
    {
        string instruction = "bash -c \"./scripts/karect_script "+file+" \"";
        cout << "Instruction: "<<instruction<<endl;
        System::execute(instruction);
        System::changeName("karect_output_file.fasta",output);
        return output;
    }
    size_t _launchCounter(string file, envi environment, bool all = false)
    {
        string instruction = "bash -c \"./scripts/dsk_script "+file+" "+to_string(KMER_SIZE)+" output.h5 output.fa "+to_string(MIN_ABUNDANCE)+" >/dev/null 2>&1\"";
        cout << "Instruction: " << instruction << endl;
        System::execute(instruction);
        cout << "Starting parse"<<endl;
        size_t total_kmers = 0;
        if (all)
            Bio::dskCounter(_kmers_map, _last_kmer, "output.fa");
        else
            total_kmers = Bio::dskParser(_envi_kmers_map, _kmers_map, "output.fa", environment);
        cout << "End Count and parsing"<<endl;
        return total_kmers;
    }
    void _checkSamples()
    {
        for (auto caso:_envi_kmers_map)
        {
            for (auto kmer_pair:caso.second)
            {
                size_t column = _kmers_map[kmer_pair.first], row = caso.first;
                _fsContext->update(column, row,
                        (double)kmer_pair.second/(double)_samples_kmers_map[caso.first], floor(caso.first / NUM_SAMPLES), 1);
            }
            cout << endl;
        }
    }
    SequenceContainer _sc;
    vector<string> _reads, _environments;
    unordered_map<envi, size_t> _samples_kmers_map;
    unordered_map<string, envi> _samples_map;
    unordered_map<Kmer, size_t> _kmers_map;
    unordered_map<envi, unordered_map<Kmer, size_t>> _envi_kmers_map;
    /*
     * Provisional
     */
    vector<Kmer> _vector_kmer_tmp;
    FSContext * _fsContext;
    ProbTable * _probTable;
    size_t _reads_processed, _last_kmer;
};
#endif