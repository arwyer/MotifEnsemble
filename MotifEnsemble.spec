/*
A KBase module: MotifEnsemble
*/

module MotifEnsemble {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_MotifEnsemble(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

    /*
    Internal workflow:
    1. Input - list of motifsets , workspace, threshold consensus
    2. Download MotifSets -> Utils function
    3. Assign motif ids by position in list

    Use refs to identify MSOs internally!
    Dictionary of motifsets key: ref, val set

    list of match sets:
    each item in the set is a tuple of (ref,index)



    for each motifset: <- enumerate to avoid duplicate
      for each motif in motifset
        for each other motifset: <- enumerate to avoid duplicate
          for each motif in other:
            compare(motif1,motif2):
              if motifs same:
                search list of sets for motif1:
                  if found add  motif2 if not in
                if not found search list of sets for motif2:
                  if found add motif1
                else add a new set with motif1 + motif2



    */

    typedef structure{
      list<string> motifset_refs;
      string workspace_name;
      float threshold;
      float proportion;
    } EnsembleParams;

    typedef structure{
      string report_name;
      string report_ref;
    } Ensemble_out;


    funcdef MotifEnsemble(EnsembleParams params)
      returns (Ensemble_out out) authentication required;
};
