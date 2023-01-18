import numpy as np
import json, os
from rdflib import Graph
from rdflib import URIRef, Literal, Namespace
from rdflib.namespace import RDF, RDFS

import urllib

import uuid
import rdflib
from SPARQLWrapper import SPARQLWrapper, JSON

from franz.openrdf.sail.allegrographserver import AllegroGraphServer
from franz.openrdf.connect import ag_connect
from franz.openrdf.rio.rdfformat import RDFFormat
from franz.openrdf.repository.repository import Repository

class Triplification_process:
    graph_proteins = Graph()
    taxon_library = {}

    def __init__(self):
        if(os.path.isfile("taxon_library.tsv")):
            f=open("taxon_library.tsv","r")
            for line in f:
                l=line.replace("\n","").split("\t")
                self.taxon_library[l[0]]=l[1]
            f.close()

        base_annotation_folder="knowledge_base_proteins"
        if(os.path.isdir(base_annotation_folder)):
            if(os.path.isfile(base_annotation_folder+'/dataset_annotation_info.ttl')):
                self.graph_proteins=Graph()
                self.graph_proteins.parse(base_annotation_folder+'/dataset_annotation_info.ttl', format="turtle")

    def find_taxon_uniprot(self, sciname):
        taxon=""
        if(not sciname.lower() in self.taxon_library.keys()):
            sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql")
            sparql.setQuery("""
                PREFIX up: <http://purl.uniprot.org/core/>
                PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
                PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                PREFIX owl: <http://www.w3.org/2002/07/owl#>
                SELECT ?taxon ?name
                WHERE
                {
                    ?taxon a up:Taxon .
                    ?taxon up:scientificName ?name .
                    filter( regex(?name, '"""+sciname+"""', 'i') ) . 
                }
            """)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()

            for result in results["results"]["bindings"]:
                _taxon=result["taxon"]["value"]
                name=result["name"]["value"]
                if(name.lower()==sciname.lower()):
                    self.taxon_library[sciname.lower()]=_taxon
                    with open("taxon_library.tsv", "a") as gf:
                        gf.write(sciname.lower()+"\t"+_taxon+"\n")
                    taxon=_taxon

        else:
            taxon= self.taxon_library[sciname.lower()]

        return taxon


    def find_protein_identifier(self, uniprotid):
        uri=""
        results = self. graph_proteins.query("""
        SELECT ?proteinid WHERE {
            ?proteinid <https://www.ypublish.info/protein_interaction_domain_ontology#hasUniprotCorrespondent> <http://purl.uniprot.org/uniprot/%s> .
        }
        """ % uniprotid)
        for row in results.result:
            uri=row[0]
            
        return uri

    def generate_info_experiment(self, config_exp, evidences_file):
        
        with open(config_exp) as json_file:  
            data = json.load(json_file)

        ppiprov = Namespace("https://www.ypublish.info/provenance_information#")
        ontoppi = Namespace("https://www.ypublish.info/protein_interaction_domain_ontology#")
        owl = Namespace("http://www.w3.org/2002/07/owl#")
        
        name_exp=data["name"].lower().replace(" ","_")

        g = Graph()

        g.add( (eval("ppiprov.exp_"+name_exp), RDF.type, ontoppi.Experiment) )
        g.add( (eval("ppiprov.exp_"+name_exp), RDFS.label, Literal(data["name"])) )
        g.add( (eval("ppiprov.exp_"+name_exp), ontoppi.hasDescription, Literal(data["description"])) )
        g.add( (eval("ppiprov.exp_"+name_exp), ontoppi.hasOwner, Literal(data["owner"])) )
        g.add( (eval("ppiprov.exp_"+name_exp), ontoppi.hasContactEmail, Literal(data["email"])) )
        #g.add( (ppiprov.exp1, ontoppi.finalPredictionMethod, Literal("PredRep Combination with Adaboost classifier")) )

        folder_files=name_exp
        if(not os.path.isdir(folder_files)):
            os.system("mkdir "+folder_files)

        f=open(folder_files+"/provenance_info.ttl","w")
        f.close()

        error_org=False
        count=1
        for d in data["datasets"]:
            name=d["name"].lower().replace(" ","_")
            g.add( (eval('ppiprov.dataset_'+name), RDF.type, ontoppi.PPIDataset) )
            g.add( (eval('ppiprov.dataset_'+name), RDFS.label, Literal(d["name"])) )
            g.add( (eval('ppiprov.dataset_'+name), ontoppi.hasDescription, Literal(d["description"])) )
            g.add( (eval('ppiprov.dataset_'+name), ontoppi.hasFolder, Literal(d["folder"])) )
            g.add( (eval('ppiprov.dataset_'+name), ontoppi.hasPrefix, Literal(d["prefix"])) )

            taxon=self.find_taxon_uniprot(d["organism"])
            if(taxon!=""):
                sciname=d["organism"].replace(" ","_").lower()
                g.add( (eval('ppiprov.organism_'+sciname), RDF.type, URIRef("http://rdf.geospecies.org/ont/geospecies#IndividualOrganism")) )
                g.add( (eval('ppiprov.organism_'+sciname), owl.sameAs, URIRef(taxon)) )
                g.add( (eval('ppiprov.dataset_'+name), ontoppi.fromOrganism, eval('ppiprov.organism_'+sciname) ) )
            else:
                error_org=True

            if(not error_org):
                g.add( (eval('ppiprov.exp_'+name_exp), ontoppi.hasDataset, eval('ppiprov.dataset_'+name)) )
                count+=1

        if(not error_org):
            with open(evidences_file) as json_file:  
                evidences = json.load(json_file)
                evidences=evidences["evidences"]
            
            count=1
            for e in evidences:
                name=e["name"].lower().replace(" ","_")
                g.add( (eval('ppiprov.evidenceMethod_'+name), RDF.type, ontoppi.EvidenceMethod) )
                g.add( (eval('ppiprov.evidenceMethod_'+name), RDFS.label, Literal(e["name"])) )
                g.add( (eval('ppiprov.evidenceMethod_'+name), ontoppi.hasDescription, Literal(e["description"])) )
                g.add( (eval('ppiprov.evidenceMethod_'+name), ontoppi.hasType, Literal(e["type_"])) )
                g.add( (eval('ppiprov.exp_'+name_exp), ontoppi.usesEvidenceMethod, eval('ppiprov.evidenceMethod_'+name)) )
                count+=1

            g.serialize(destination=folder_files+'/provenance_info.ttl', format='turtle')
        else:
            print("Error: The organism name given is wrong and does not match with any known organism")

        return error_org

    def generate_annotation(self, config_exp):
        base_annotation_folder="knowledge_base_proteins"
        if(not os.path.isdir(base_annotation_folder)):
            os.system("mkdir "+base_annotation_folder)

        ontoppi = Namespace("https://www.ypublish.info/protein_interaction_domain_ontology#")

        with open(config_exp) as json_file:  
            data = json.load(json_file)

        predicates_features=['ontoppi.hasGO_cc_annotation','ontoppi.hasGO_mf_annotation','ontoppi.hasGO_bp_annotation','ontoppi.hasKo_annotation','ontoppi.hasPfam_annotation']
        for d in data["datasets"]:
            name=d["name"].lower().replace(" ","_")
            
            cnt=0
            
            if( os.path.isfile(d["folder"]+"new_info_proteins_false.txt") and os.path.isfile(d["folder"]+"new_info_proteins_positive.txt") ):
                proteins=[]
                ds=['false','positive']
                for d_ in ds:
                    f=open(d["folder"]+"new_info_proteins_"+d_+".txt","r")
                    for line in f:
                        uniqueid=str(uuid.uuid4())

                        l=line.replace("\n","").split("\t")
                        if(not l[0] in proteins):
                            proteins.append(l[0])
                            
                            id=self.find_protein_identifier(l[0])
                            if(id==""):
                                self. graph_proteins.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+uniqueid), RDF.type, ontoppi.PairComponent) )
                                self. graph_proteins.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+uniqueid), ontoppi.hasUniprotCorrespondent, URIRef("http://purl.uniprot.org/uniprot/"+l[0])) )
                                
                                c=0
                                for info in l[1:-1]:
                                    if(info!="None" and info!=""):
                                        data=info.split(" ")
                                        for dt in data:
                                            if(dt!="None" and dt!=""):
                                                self. graph_proteins.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+uniqueid), eval(predicates_features[c]), Literal(dt)) )
                                    c+=1

                            #if(cnt%5000==0 and cnt!=0):
                            #    g.serialize(destination=base_annotation_folder+'/dataset_part-'+str(cnt)+'_'+str(i)+'_annotation_info.ttl', format='turtle')
                            #    g=Graph()
                            cnt+=1

                    f.close()
                self. graph_proteins.serialize(destination=base_annotation_folder+'/dataset_annotation_info.ttl', format='turtle')
            #else:
            #    print("Error the files with protein functional features were not found in datasets folders")
            
    def generate_results(self, config_exp, evidences_file):
        
        with open(config_exp) as json_file:  
            data = json.load(json_file)

        name_exp=data["name"].lower().replace(" ","_")

        g = Graph()

        folder_files=name_exp+"/"
        if(not os.path.isdir(folder_files)):
            os.system("mkdir "+folder_files)

        ontoppi = Namespace("https://www.ypublish.info/protein_interaction_domain_ontology#")
        ppiprov = Namespace("https://www.ypublish.info/provenance_information#")
        biopax = Namespace("http://www.biopax.org/release/biopax-level2.owl#")
        
        with open(evidences_file) as json_file:  
            evidences = json.load(json_file)
            evidences=evidences["evidences"]

        literature_support=False
        for e in evidences:
            if("literature_support" in e.keys()):
                if(e["literature_support"]):
                    literature_support=True
        
        pairs={}
        for d in data["datasets"]:
            name=d["name"].lower().replace(" ","_")
            pairs[name]=[]
            f=open(d["folder"]+d['prefix']+"dataset.txt","r")
            for line in f:
                l=line.replace("\n","").split("\t")
                pairs[name].append([l[0],l[1]])
            f.close()

        for d in data["datasets"]:
            if(literature_support and os.path.isfile(d["folder"]+"literature_link.tsv") ):
                pubmedids={}
                f=open(d["folder"]+"literature_link.tsv", "r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    pubmedids[l[0]+"-"+l[1]]=l[2].split(" ")
                f.close()

            name=d["name"].lower().replace(" ","_")
            f=open(folder_files+'dataset_'+name+'_results_info.ttl',"w")
            f.close()

            f=open(folder_files+'interactions_'+name+'.tsv',"w")
            f.close()
            
            g=Graph()
            c=0
            f=open(d["folder"]+"complete_dataset_ppi.tsv","r")
            for line in f:
                uniqueid=str(uuid.uuid4())
                l=line.replace("\n","").split("\t")

                p1=pairs[name][c][0]
                p2=pairs[name][c][1]
                id1=self.find_protein_identifier(p1)
                id2=self.find_protein_identifier(p2)
                g.add( (URIRef('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid), RDF.type, biopax.Interaction ) )

                if(id1==""):
                    id_protein=str(uuid.uuid4())
                    self.graph_proteins.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+id_protein), RDF.type, ontoppi.PairComponent) )
                    g.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+id_protein), ontoppi.hasUniprotCorrespondent, URIRef("http://purl.uniprot.org/uniprot/"+p1)) )
                    id1='https://www.ypublish.info/protein_annotation_information#protein_'+id_protein
                
                if(id2==""):
                    id_protein=str(uuid.uuid4())
                    self.graph_proteins.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+id_protein), RDF.type, ontoppi.PairComponent) )
                    g.add( (URIRef('https://www.ypublish.info/protein_annotation_information#protein_'+id_protein), ontoppi.hasUniprotCorrespondent, URIRef("http://purl.uniprot.org/uniprot/"+p2)) )
                    id2='https://www.ypublish.info/protein_annotation_information#protein_'+id_protein
                    
                g.add( (URIRef('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid), ontoppi.participant1, URIRef(id1) ) )
                g.add( (URIRef('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid), ontoppi.participant2, URIRef(id2) ) )
                
                g.add( (URIRef('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid), ontoppi.belongsTo, eval('ppiprov.dataset_'+name) ) )

                with open(folder_files+'interactions_'+name+'.tsv', 'a') as gf:
                    gf.write('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid+"\t"+p1+"\t"+p2+"\n")
                
                cc=0
                data=[]
                for info in l:
                    evname=evidences[cc]["name"].lower().replace(" ","_")
                    g.add( (URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid), RDF.type, ontoppi.EvidenceDecision ) )
                    g.add( (URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid), ontoppi.basedOn, eval('ppiprov.evidenceMethod_'+evname) ) )
                    g.add( (URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid), ontoppi.predictedValue, Literal(float(info)) ) )
                    
                    if("literature_support" in evidences[cc].keys()):
                        if(evidences[cc]["literature_support"] and os.path.isfile(d["folder"]+"literature_link.tsv")):
                            if( (p1+"-"+p2) in pubmedids.keys()):
                                for pubid in pubmedids[p1+"-"+p2]:
                                    g.add( (URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid), ontoppi.supportingArticle, Literal(pubid) ) )
                            else:
                                if( (p2+"-"+p1) in pubmedids.keys()):
                                    for pubid in pubmedids[p2+"-"+p1]:
                                        g.add( (URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid), ontoppi.supportingArticle, Literal(pubid) ) )
                    
                    g.add( (URIRef('https://www.ypublish.info/prediction_results_information#interaction_'+uniqueid), ontoppi.hasScore, URIRef('https://www.ypublish.info/prediction_results_information#evidenceDecision_dataset_'+name+"_"+evname+"_int"+uniqueid) ) )
                    cc+=1

                #if(c%2000==0 and c!=0): 
                #    g.serialize(destination=folder_files+'dataset_'+name+'_results_info.ttl', format='turtle')
                #    g=Graph()
                c+=1
            f.close()
            g.serialize(destination=folder_files+'dataset_'+name+'_results_info.ttl', format='turtle')

        base_annotation_folder="knowledge_base_proteins"
        if(not os.path.isdir(base_annotation_folder)):
            os.system("mkdir "+base_annotation_folder)
        self.graph_proteins.serialize(destination=base_annotation_folder+'/dataset_annotation_info.ttl', format='turtle')
    
    def intersection(self, lst1, lst2): 
        return [item for item in lst1 if item in lst2]
    
    def execute_data_fusion(self, config_evidence_file):
        owl = Namespace("http://www.w3.org/2002/07/owl#")

        information=[]
        f=open(config_evidence_file,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            information.append([l[0], l[1]])
        f.close()

        if(not os.path.isdir("mappings")):
            os.system("mkdir mappings")

        organisms=[]
        for i in information:
            with open(i[0]) as json_file:  
                data = json.load(json_file)

            for d in data["datasets"]:
                sciname=d["organism"].lower().replace(" ","_")
                if(not sciname in organisms):
                    organisms.append(sciname)
            
        for o in organisms:
            g=Graph()
            pairs_datasets={}
            identifiers={}
            for i in information:
                with open(i[0]) as json_file:  
                    data = json.load(json_file)

                name_exp=data["name"].lower().replace(" ","_")
                folder_files=name_exp+"/"

                for d in data["datasets"]:
                    sciname=d["organism"].lower().replace(" ","_")
                    name=d["name"].lower().replace(" ","_")
                    
                    if(sciname==o):
                        pairs_datasets[name]=[]
                        identifiers[name]={}
                        f=open(folder_files+'interactions_'+name+'.tsv', 'r')
                        for line in f:
                            l=line.replace("\n","").split("\t")
                            sort=[l[1], l[2]]
                            sort.sort()
                            pairs_datasets[name].append(sort)
                            identifiers[name][sort[0]+"-"+sort[1]]=l[0]
                        f.close()
            passed=[]
            for ds in pairs_datasets.keys():
                for ds2 in pairs_datasets.keys():
                    if(ds!=ds2 and not(ds+"-"+ds2 in passed)):
                        passed.append(ds+"-"+ds2)
                        passed.append(ds2+"-"+ds)
                        intersection = self.intersection(pairs_datasets[ds], pairs_datasets[ds2])
                        for inter in intersection:
                            id_=inter[0]+"-"+inter[1]
                            g.add( (URIRef(identifiers[ds][id_]), owl.sameAs, URIRef(identifiers[ds2][id_])) )
                         

            g.serialize(destination='mappings/mappings_'+o+'.ttl', format='turtle')

    def export_to_allegrograph(self, config_evidence_file):
        try:
            print("Getting information about allegrograph user")
            AGRAPH_HOST = os.environ.get('AGRAPH_HOST')
            AGRAPH_PORT = int(os.environ.get('AGRAPH_PORT', '10035'))
            AGRAPH_USER = os.environ.get('AGRAPH_USER')
            AGRAPH_PASSWORD = os.environ.get('AGRAPH_PASSWORD')
            
            print("Initializing repository")
            server = AllegroGraphServer(AGRAPH_HOST, AGRAPH_PORT, AGRAPH_USER, AGRAPH_PASSWORD)
            catalog = server.openCatalog('')

            mode = Repository.RENEW
            repo = catalog.getRepository('ppitriplificator', mode)
            repo.initialize()
            conn = repo.getConnection()

            try:
                information=[]
                f=open(config_evidence_file,"r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    information.append([l[0], l[1]])
                f.close()

                organisms=[]
                for i in information:
                    with open(i[0]) as json_file:  
                        data = json.load(json_file)

                    name_exp=data["name"].lower().replace(" ","_")
                    folder_files=name_exp

                    for d in data["datasets"]:
                        sciname=d["organism"].lower().replace(" ","_")
                        if(not sciname in organisms):
                            organisms.append(sciname)

                print("Loading ontology ontoppi for reasoning inference")
                conn.addFile('ontoppi.ttl', None, format=RDFFormat.TURTLE, context=None)
    
                print("Loading knowledge base about proteins")
                base_annotation_folder="knowledge_base_proteins"
                conn.addFile(base_annotation_folder+'/dataset_annotation_info.ttl', None, format=RDFFormat.TURTLE, context=None)

                print("Loading provenance information about the experiments and score results for the interactions")
                for i in information:
                    with open(i[0]) as json_file:  
                        data = json.load(json_file)
                    
                    name_exp=data["name"].lower().replace(" ","_")
                    folder_files=name_exp+"/"
                    if(os.path.isdir(folder_files)):
                        conn.addFile(folder_files+"provenance_info.ttl", None, format=RDFFormat.TURTLE, context=None)

                        for d in data["datasets"]:
                            name=d["name"].lower().replace(" ","_")
                            conn.addFile(folder_files+'dataset_'+name+'_results_info.ttl', None, format=RDFFormat.TURTLE, context=None)

                print("Loading interaction mappings between different datasets for each organism")
                for o in organisms:
                    conn.addFile('mappings/mappings_'+o+'.ttl', None, format=RDFFormat.TURTLE, context=None)
            except:
                print("Error: The files are not in the correct format or you did not run all steps required to execute this task")
                
        except:
            print("Error when trying to communicate with allegrograph server\n\
                Visit this pages to install it:\n\
                \tServer: https://franz.com/agraph/support/documentation/current/server-installation.html\
                \tClient: https://franz.com/agraph/support/documentation/current/python/install.html\
                \tStart server: ~/allegrograph/bin/agraph-control --config ~/allegrograph/lib/agraph.cfg start\
                \tDefine the folllowing environment variables in ~/.bashrc : AGRAPH_HOST, AGRAPH_PORT, AGRAPH_USER and AGRAPH_PASSWORD")



class Running_config:
    a=Triplification_process()

    def run_step1(self, config_exp, evidences_file):
        print("Running step 1")
        flag=self.a.generate_info_experiment(config_exp, evidences_file)
        return flag

    def run_step2(self, config_exp):
        print("Running step 2")
        self.a.generate_annotation(config_exp)

    def run_step3(self, config_exp, evidences_file):
        print("Running step 3")
        self.a=Triplification_process()
        self.a.generate_results(config_exp, evidences_file)

    def run_step4(self, config_evidence_file):
        print("Running step 4")
        self.a.execute_data_fusion(config_evidence_file)

    def run_step5(self, config_evidence_file):
        information=[]
        f=open(config_evidence_file,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            information.append(l)
        f.close()

        print("Running all the processes")
        for i in information:
            print("\tRunning now configuration for ", i[0])
            self.run_step1(i[0], i[1])
            self.run_step2(i[0])
            self.run_step3(i[0], i[1])

        self.run_step4(config_evidence_file)

    def run_step6(self, config_evidence_file):
        print("Running Exporting to allegrograph")
        self.a.export_to_allegrograph(config_evidence_file)

    def run(self, args):
        error=False
        run=0
        if(args.running_type=="" ):
            run=0
        else:
            if(args.running_type in [0,1,2,3,4,5,6]):
                run=args.running_type
            else:
                print("Error: invalid choice")
                error=True

        if(not error):
            if(run==4 or run==5 or run==6):
                if(args.file_config_evidence!=""):
                    if(run==4):
                        self.run_step4(args.file_config_evidence)
                    if(run==5):
                        self.run_step5(args.file_config_evidence)
                    if(run==6):
                        self.run_step6(args.file_config_evidence)
                else:
                    print("Error: You have to specify a valid file in tsv with the addresses of configuration and evidences files by line")
            else:
                if(args.file_experiment_config!="" and os.path.isfile(args.file_experiment_config)):
                    if(run==0 or run==1 or run==3):
                        if(args.file_config_evidence==""):
                            print("Error: you have to give the file with the evidences specification")
                        else:
                            
                            if(run==0):
                                print("Running step 1")
                                flag=self.run_step1(args.file_experiment_config , args.file_evidence_info)

                                if(not flag):
                                    print("Running step 2")
                                    self.run_step2(args.file_experiment_config)

                                    print("Running step 3")
                                    self.run_step3(args.file_experiment_config, args.file_evidence_info)

                            if(run==1):
                                self.run_step1(args.file_experiment_config, args.file_evidence_info)

                            if(run==3):
                                self.run_step3(args.file_experiment_config, args.file_evidence_info)

                    else:
                        self.run_step2(args.file_experiment_config)
                else:
                    print("Error: You have to specify a valid experiment configuration file")

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' PPIDataTriplificator - Generating protein interactions information for the linked data cloud', formatter_class=RawTextHelpFormatter)
parser.add_argument("-rt", "--running_type", action="store", help="0 - Generate the descriptions for all the protein interaction steps of an experiment  (run steps 1, 2 and 3)\n\
    1 - Generate triples just about data provenance \n\
    2 - Generate triples just for protein functional annotations\n\
    3 - Generate triples just for the score results of each evidence\n\
    4 - Execute data fusion\n\
    5 - Generate descriptions and execute data fusion (run steps 1, 2, 3 and 4)\n\
    6 - Export to allegrograph server", type=int)
parser.add_argument("-fec", "--file_experiment_config", action="store", help="File with the experiment configuration in json format")
parser.add_argument("-fev", "--file_evidence_info", action="store", help="File with the evidence methods used and their information in json format")
parser.add_argument("-fcv", "--file_config_evidence", action="store", help="File with the experiment and evidence methods files addresses in tsv format")
args = parser.parse_args()
r=Running_config()
r.run(args)
