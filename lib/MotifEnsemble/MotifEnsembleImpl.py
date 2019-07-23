# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from installed_clients.KBaseReportClient import KBaseReport
from MotifEnsemble.Utils.DownloadMotifSets import DownloadMotifSets
from MotifEnsemble.Utils.FastaUtils import FastaUtils
from installed_clients.DataFileUtilClient import DataFileUtil
from copy import deepcopy
import uuid
from MotifEnsemble.Utils.MakeNewReport import MakeNewReport
#END_HEADER


class MotifEnsemble:
    '''
    Module Name:
    MotifEnsemble

    Module Description:
    A KBase module: MotifEnsemble
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbasecollaborations/MotifEnsemble.git"
    GIT_COMMIT_HASH = "d414d0eee0a62d4c12e6368ce14571900c291b53"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def run_MotifEnsemble(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_MotifEnsemble
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_MotifEnsemble

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_MotifEnsemble return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def MotifEnsemble(self, ctx, params):
        """
        :param params: instance of type "EnsembleParams" (Internal workflow:
           1. Input - list of motifsets , workspace, threshold consensus 2.
           Download MotifSets -> Utils function 3. Assign motif ids by
           position in list Use refs to identify MSOs internally! Dictionary
           of motifsets key: ref, val set list of match sets: each item in
           the set is a tuple of (ref,index) for each motifset: <- enumerate
           to avoid duplicate for each motif in motifset for each other
           motifset: <- enumerate to avoid duplicate for each motif in other:
           compare(motif1,motif2): if motifs same: search list of sets for
           motif1: if found add  motif2 if not in if not found search list of
           sets for motif2: if found add motif1 else add a new set with
           motif1 + motif2) -> structure: parameter "motifset_refs" of list
           of String, parameter "workspace_name" of String, parameter
           "threshold" of Double, parameter "proportion" of Double
        :returns: instance of type "Ensemble_out" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: out
        #BEGIN MotifEnsemble
        #TODO: ERROR CHECK (MULTIPLE MOTIFSETS, NONEMPTY, SSREF are the same, etc.)

        dms=DownloadMotifSets()
        MotifSetDict = dms.DownloadMotifSet(params['motifset_refs'],self.callback_url)

        matchSets = []
        threshold = float(params['threshold'])
        fmu=FastaUtils()
        for i,MSR1 in enumerate(MotifSetDict.keys()):
            for j,motif1 in enumerate(MotifSetDict[MSR1]['Motifs']):
                for k,MSR2 in enumerate(MotifSetDict.keys()):
                    if k > i:
                        for l,motif2 in enumerate(MotifSetDict[MSR2]['Motifs']):
                            print(motif1)
                            print(motif2)
                            print(threshold)
                            if fmu.CompareMotifsBP(motif1,motif2,threshold):
                                found1 = False
                                found2 = False
                                index1 = -1
                                index2 = -1
                                for m,mset in enumerate(matchSets):
                                    if (MSR1,j) in mset:
                                        found1 = True
                                        index1 = m
                                    if(MSR2,l) in mset:
                                        found2 = True
                                        index2 = m
                                if not found1 and found2:
                                    matchSets[index2].add((MSR1,j))
                                elif not found2 and found1:
                                    matchSets[index1].add((MSR2,l))
                                elif found1 and found2:
                                    if index1 != index2:
                                        matchSets[index1].union(matchSets[index2])
                                        matchSets.pop(index2)
                                else:
                                    matchSets.append(set([(MSR1,j),(MSR2,l)]))
        numMotifSets = len(params['motifset_refs'])
        threshold = float(params['proportion'])
        KeepSets = []
        print('NUM MATCHSETS********')
        print(len(matchSets))
        for i,mset in enumerate(matchSets):
            uniqueRefs = {}
            for tuple in mset:
                if tuple[0] not in uniqueRefs:
                    uniqueRefs[tuple[0]] = tuple[0]
            if float(len(uniqueRefs.keys()))/numMotifSets >= threshold:
                KeepSets.append(i)
        print(len(KeepSets))

        ESO = {}
        for ref in MotifSetDict:
            ESO['Condition'] = MotifSetDict[ref]['Condition']
            ESO['SequenceSet_ref'] = MotifSetDict[ref]['SequenceSet_ref']
            ESO['Alphabet'] = deepcopy(MotifSetDict[ref]['Alphabet'])
            ESO['Background'] = deepcopy(MotifSetDict[ref]['Background'])
            break
        ESO['Motifs'] = []
        #Add motifs
        for keep in KeepSets:
            motif = fmu.merge(matchSets[keep],MotifSetDict)
            ESO['Motifs'].append(deepcopy(motif))


        #upload new MSO
        dfu = DataFileUtil(self.callback_url)
        save_objects_params = {}
        save_objects_params['id'] = dfu.ws_name_to_id(params['workspace_name'])
        save_objects_params['objects'] = [{'type': 'KBaseGeneRegulation.MotifSet' , 'data' : ESO , 'name' : 'EnsembleMotifSet'}]

        info = dfu.save_objects(save_objects_params)[0]
        obj_ref = "%s/%s/%s" % (info[6], info[0], info[4])
        htmlDir = self.shared_folder + '/ensemble_html'
        os.mkdir(htmlDir)
        mr=MakeNewReport()
        mr.MakeReport(htmlDir,ESO)


        try:
            html_upload_ret = dfu.file_to_shock({'file_path': htmlDir ,'make_handle': 0, 'pack': 'zip'})
        except:
            raise ValueError ('error uploading HTML file to shock')


        reportName = 'MEMEMotifFinder_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [{'ref' : obj_ref, 'description' : 'Motif Set generated by MEME'}],
                     'message': '',
                     'direct_html': None,
                     'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }


        # attach to report obj
        reportObj['direct_html'] = ''
        reportObj['direct_html_link_index'] = 0
        reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                    'name': 'index.html',
                                    'label': 'Save promoter_download.zip'
                                    }
                                   ]


        report = KBaseReport(self.callback_url, token=ctx['token'])
        report_info = report.create_extended_report(reportObj)
        out = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END MotifEnsemble

        # At some point might do deeper type checking...
        if not isinstance(out, dict):
            raise ValueError('Method MotifEnsemble return value ' +
                             'out is not type dict as required.')
        # return the results
        return [out]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
