import json
import os
import numpy as np

class PPIPreTriplification:
    def get_information_dip(self, folder):
        if(os.path.isfile(folder+"dip_info.txt")):
            f=open(folder+"literature_link.tsv","w")
            f.close()

            i=0
            f=open(folder+"dip_info.txt","r")
            for line in f:
                l=line.split("\t")
                if(i>0):
                    p1=""
                    p1l=l[0].split("|")
                    for p in p1l:
                        if(p.find("uniprotkb")!=-1):
                            p1=p.replace("uniprotkb:","")
                    p2=""
                    p2l=l[1].split("|")
                    for p in p2l:
                        if(p.find("uniprotkb")!=-1):
                            p2=p.replace("uniprotkb:","")

                    exp_methods=l[6].split("|")
                    mi_ids=[]
                    names=[]
                    for e in exp_methods:
                        aux=e.replace(")","").split("(")
                        mi_ids.append(aux[0])
                        names.append(aux[1]) 

                    papers=l[8].split("|")
                    pmids=[]
                    for p in papers:
                        aux=p.replace("pubmed:","")
                        if(aux.find("DIP-")==-1 and not(aux in pmids)):
                            pmids.append(aux)

                    with open(folder+"literature_link.tsv","a")as gf:
                        gf.write("%s\t%s\t%s\n" %(p1, p2, (" ".join(pmids)) ) )
                        #gf.write(p1+";"+p2+";"+("|".join(mi_ids))+";"+("|".join(names))+";"+("|".join(pmids))+"\n")

                i+=1
            f.close()

    def string_case(self, config_exp, organism=None):
        with open(config_exp) as json_file:  
            data = json.load(json_file)
            
        condition=True
        for d in data["datasets"]:
            if(organism!=None):
                condition=(organism==d["organism"])

            if(condition):
                f=open(d["folder"]+"complete_dataset_ppi.tsv","w")
                f.close()

                f=open(d["folder"]+"string_dataset.tsv","w")
                f.close()

                mapp={}
                f=open(d["folder"]+"string_mapping.tsv","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    mapp[l[1].split("|")[0]]=l[2]
                f.close()
                
                pairs=[]
                f=open(d["folder"]+"predicted_as_positive.tsv","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    pairs.append([l[0],l[1]])
                f.close()

                for p in pairs:
                    if(p[0] in mapp.keys() and p[1] in mapp.keys()):
                        #p1=mapp[p[0]]
                        #p2=mapp[p[1]]

                        p1=p[0]
                        p2=p[1]
                        
                        with open(d["folder"]+"string_dataset.txt","a") as g:
                            g.write("%s\t%s\t1\n" %(p[0], p[1]) )

                        os.system("grep "+p1+" "+d["folder"]+"string_unique_interactions.tsv | grep "+p2+" > "+d["folder"]+"filter.txt")
                        f=open(d["folder"]+"filter.txt","r")
                        for line in f:
                            r=line.replace("\n","").split("\t")
                            vector=r[2:]
                            l="\t".join(vector)
                            with open(d["folder"]+"complete_dataset_ppi.tsv","a") as g:
                                g.write(l+"\n")

                            break
                        f.close()

    def hint_case(self, config_exp, organism=None):
        with open(config_exp) as json_file:  
            data = json.load(json_file)
            
        condition=True
        for d in data["datasets"]:
            if(organism!=None):
                condition=(organism==d["organism"])

            if(condition):
                f=open(d["folder"]+"complete_dataset_ppi.tsv","w")
                f.close()

                f=open(d["folder"]+"hint_dataset.txt","w")
                f.close()

                c=0
                f=open(d["folder"]+"hint_unique_interactions.txt","r")
                for line in f:
                    l=line.split("\t")
                    with open(d["folder"]+"complete_dataset_ppi.tsv","a") as g:
                        g.write("1\n")
                    
                    with open(d["folder"]+"hint_dataset.txt","a") as g:
                        g.write("%s\t%s\t1\n" %(l[0], l[1]) )

                    c+=1
                f.close()

            self.get_information_dip(d["folder"])

    def predppi_case(self, config_exp, organism=None):
        with open(config_exp) as json_file:  
            data = json.load(json_file)
            
        condition=True
        for d in data["datasets"]:
            if(organism!=None):
                condition=(organism==d["organism"])

            if(condition):
                f=open(d["folder"]+"complete_dataset_ppi.tsv","w")
                f.close()

                f=open(d["folder"]+"predicted_as_positive.tsv","w")
                f.close()
                pairs=[]
                f=open(d["folder"]+d['prefix']+"dataset.txt","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    pairs.append([l[0],l[1]])
                f.close()

                final_score=np.load(d["folder"]+"final_score.npy")
                c=0
                f=open(d["folder"]+"dataset_ppi.txt","r")
                for line in f:
                    l=line.split(" ")
                    l=l[1:-1]
                    l.append(str(final_score[c]))
                    l="\t".join(l)
                    with open(d["folder"]+"complete_dataset_ppi.tsv","a") as g:
                        g.write(l+"\n")
                    
                    if(final_score[c]>0.8):
                        with open(d["folder"]+"predicted_as_positive.tsv","a") as g:
                            g.write( ("\t".join(pairs[c]))+"\n")
                    c+=1
                f.close()
            
class Running_config:

    def run_step1(self, config_exp, organism):
        a=PPIPreTriplification()
        a.predppi_case(config_exp, organism)

    def run_step2(self, config_exp, organism):
        a=PPIPreTriplification()
        a.string_case(config_exp, organism)

    def run_step3(self, config_exp, organism):
        a=PPIPreTriplification()
        a.hint_case(config_exp, organism)

    def run(self, args):
        run=0
        if(args.running_type=="" ):
            run=0
        else:
            if(args.running_type in [1,2,3]):
                run=args.running_type
            else:
                print("Error: invalid choice")

        if(args.file_experiment_config!="" and os.path.isfile(args.file_experiment_config)):
            if(run==1):
                self.run_step1(args.file_experiment_config, args.organism)

            if(run==2):
                self.run_step2(args.file_experiment_config, args.organism)

            if(run==3):
                self.run_step3(args.file_experiment_config, args.organism)
        else:
            print("Error: You have to specify a valid experiment configuration file")

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' PPIPreTriplification - Tool to organize data from HINT, PredPPI and String to triplify', formatter_class=RawTextHelpFormatter)
parser.add_argument("-rt", "--running_type", action="store", help="\
    1 - Prepare data for PredPPI \n\
    2 - Prepare data for String\n\
    3 - Prepare data for HINT", type=int)
parser.add_argument("-fec", "--file_experiment_config", action="store", help="File with the experiment configuration in json format")
parser.add_argument("-org", "--organism", action="store", help="Organism of interest")
args = parser.parse_args()
r=Running_config()
r.run(args)

