import os
import requests
import mygene
import random
import json
import urllib.request

class HostPathogenBacteria:
    def get_from_hgnc(self, gene, org):
        r=requests.get("https://www.uniprot.org/uniprot/?query=gene:"+gene+"+AND+organism:"+org+"&format=tab")
        for line in r.text.strip().split("\n"):
            l = line.strip().split("\t")
            if(l[0]!="Entry"):
                return l[0]
        return ''
    
    def get_from_refseq(self, gene):
        mg = mygene.MyGeneInfo()
        r=mg.querymany([gene], scopes='refseq', fields='symbol,name,uniprot')
        if('uniprot' in r[0]):
            return r[0]['uniprot']['Swiss-Prot']
        return ''
        
    def build_mapping_training_data(self):
        ok=set()
        f=open("host_pathogen_interaction/training_dataset/training_positive.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            ide=l[0]+"-"+l[1]
            if(not ide in ok):
                ok.add(ide)
        f.close()
        
        mapp={}
        for f in os.listdir("host_pathogen_interaction/"):
            if(f.endswith("csv")):
                x=open("host_pathogen_interaction/"+f,"r")
                for line in x:
                    l=line.replace("\n","").split("\t")
                    if(l[0].find("_")==-1):
                        p1=self.get_from_hgnc(l[0], "83333")
                    else:
                        p1=self.get_from_refseq(l[0])
                        
                    if(l[1].find("_")==-1):
                        p2=self.get_from_hgnc(l[1], "9606")
                    else:
                        p2=self.get_from_refseq(l[1])
                        
                    if(p1!='' and p2!=''):
                        if(not p1 in mapp.keys()):
                            mapp[p1]=[l[0], 'ecoli']
                        if(not p2 in mapp.keys()):
                            mapp[p2]=[l[1], 'human']
                            
                        ide=p1+"-"+p2
                        if(not ide in ok):
                            ok.add(ide)
                            with open("host_pathogen_interaction/training_dataset/training_positive.txt", "a") as g:
                                g.write("%s\t%s\t1\n" %(p1, p2) )
                        
                x.close()
        
        h=open("mapping_genes.tsv","w")
        h.write("gene\tuniprot\torganism\n")
        for k in mapp:
            h.write("%s\t%s\t%s\n" %(k, mapp[k][0], mapp[k][1]) )
        h.close()
        
    def complement_positive(self):
        g=open("host_pathogen_interaction/training_positive_hpidb.txt", "w")
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for f in folders:
            i=0
            f=open(f+"/ppis.txt", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(i>0):
                    id1=""
                    id2=""
                    p1=l[0].split("|")
                    for id in p1:
                        if(id.startswith("uniprotkb")):
                            id1=id.split(":")[1]
                            
                    p2=l[1].split("|")
                    for id in p2:
                        if(id.startswith("uniprotkb")):
                            id2=id.split(":")[1]
                            
                    if(id1!="" and id2!=""):
                        with open("host_pathogen_interaction/training_dataset/training_positive.txt", "a") as fg:
                            fg.write("%s\t%s\t1\n" %(id1, id2) )
                        g.write("%s\t%s\t1\n" %(id1, id2) )
                i+=1
            f.close()
        g.close()
    
    def build_train_individual_datasets(self):
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for fo in folders:
            if(os.path.isdir("host_pathogen_interaction/"+fo+"_ds/")):
                os.system("rm -rf host_pathogen_interaction/"+fo+"_ds")
            os.system("mkdir host_pathogen_interaction/"+fo+"_ds")
            
            params={}
            params['name']="Experiment for context-based PPI network construction"
            params['description']="experiment train "+fo
            params['email']="ycfrenchgirl2@gmail.com"
            params['owner']="Yasmmin"
            ds={}
            ds['name'] = "predprin_"+fo
            ds["prefix"]="training_"
            ds["folder"]="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/host_pathogen_interaction/"+fo+"_ds/"
            ds["description"]="training host pathogen dataset for models, labels were randomly selected"
            params["datasets"]=[ds]
            with open('predprin/params_train_hostpathogen_'+fo+'.json', 'w') as fp:
                json.dump(params, fp)
                
            np=0
            pos=set()
            # positive dataset
            i=0
            g=open("host_pathogen_interaction/"+fo+"_ds/training_dataset.txt", "w")
            f=open(fo+"/ppis.txt", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(i>0):
                    id1=""
                    id2=""
                    p1=l[0].split("|")
                    for id in p1:
                        if(id.startswith("uniprotkb")):
                            id1=id.split(":")[1].split("-")[0]
                    
                    if(id1=="" and len(l)>17):
                        if(l[15].startswith("UNIPROT_AC")):
                            id1=l[15].split(":")[1]
                                          
                    p2=l[1].split("|")
                    for id in p2:
                        if(id.startswith("uniprotkb")):
                            id2=id.split(":")[1].split("-")[0]
                            
                    if(id2=="" and len(l)>17):
                        if(l[16].startswith("UNIPROT_AC")):
                            id2=l[16].split(":")[1]   
                            
                    if(id1!="" and id2!=""):
                        pos.add(id1+"-"+id2)
                        pos.add(id2+"-"+id1)
                        
                        with open("host_pathogen_interaction/"+fo+"_ds/training_positive.txt", "a") as fg:
                            fg.write("%s\t%s\t1\n" %(id1, id2) )
                        g.write("%s\t%s\t1\n" %(id1, id2) )
                        np+=1
                i+=1
            f.close()
            print('-----'+fo, np)
            
            # false dataset
            human=[]
            f=open("host_pathogen_interaction/hsapiens_random_ds_200000_false.txt", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(not l[0] in human):
                    human.append(l[0])
                if(not l[1] in human):
                    human.append(l[1])
            f.close()
            
            nf=0
            pat=[]
            if(fo=='ecoli'):
                f=open("host_pathogen_interaction/ecoli_random_ds_200000_false.txt", "r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    if(not l[0] in pat):
                        pat.append(l[0])
                    if(not l[1] in pat):
                        pat.append(l[1])
                f.close()
            else:
                i=0
                f=open(fo+"/proteins.tab", "r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    if(i>0):
                        if(not l[0] in pat):
                            pat.append(l[0])
                    i+=1
                f.close()
            
            new_=[]    
            for i in range(np):
                p1=random.choice(human)
                p2=random.choice(pat)
                        
                while (p1+"-"+p2 in pos and [p1, p2] in new_):
                    p2=random.choice(pat)
                    
                new_.append( [p1, p2] )    
                with open("host_pathogen_interaction/"+fo+"_ds/training_false.txt", "a") as fg:
                    fg.write("%s\t%s\t0\n" %(p1, p2) )    
                g.write("%s\t%s\t0\n" %(p1, p2) )
            g.close()
            
    def run_train_ppi(self):
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        folders=['staphylococcus_aureus','pseudomonas_aeruginosa']
        folders=['staphylococcus_aureus']
        for fo in folders:
            f=open("execute_predppi.sh","w")
            f.write("#!/bin/bash\n")
            f.write("cd predprin/\n")
            f.write("python3 -m luigi --module main RunPPIExperiment --parameters-file params_train_hostpathogen_"+fo+".json --mode train --model None --workers 1")
            f.close()

            os.system("rm predprin/run_experiment.txt")
            print("Executing ppi prediction with host pathogen pairs "+fo)
            os.system('bash execute_predppi.sh')
            
    def build_false_dataset(self):
        pos=set()
        total_positive=0
        with open("host_pathogen_interaction/training_dataset/training_positive.txt", 'r') as fp:
            for line in fp:
                l=line.replace("\n","").split("\t")
                if(not l[0]+"-"+l[1] in pos):
                    pos.add(l[0]+"-"+l[1])
                    pos.add(l[1]+"-"+l[0])
                total_positive+=1
    
        human=[]
        f=open("host_pathogen_interaction/hsapiens_random_ds_200000_false.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(not l[0] in human):
                human.append(l[0])
            if(not l[1] in human):
                human.append(l[1])
        f.close()
        
        ecoli=[]
        f=open("host_pathogen_interaction/ecoli_random_ds_200000_false.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(not l[0] in ecoli):
                ecoli.append(l[0])
            if(not l[1] in ecoli):
                ecoli.append(l[1])
        f.close()
        
        pseudo=[]
        i=0
        f=open("pseudomonas_aeruginosa/proteins.tab", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(i>0):
                if(not l[0] in pseudo):
                    pseudo.append(l[0])
            i+=1
        f.close()
        
        myco=[]
        i=0
        f=open("mycobacterium_tuberculosis/proteins.tab", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(i>0):
                if(not l[0] in myco):
                    myco.append(l[0])
            i+=1
        f.close()
        
        stap=[]
        i=0
        f=open("staphylococcus_aureus/proteins.tab", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(i>0):
                if(not l[0] in stap):
                    stap.append(l[0])
            i+=1
        f.close()
        
        j=1
        g=open("host_pathogen_interaction/training_dataset/training_false.txt", "w")
        for i in range(total_positive):
            p1=random.choice(human)
            
            if(j==1):
                p2=random.choice(ecoli)
            if(j==2):
                p2=random.choice(myco)
            if(j==3):
                p2=random.choice(pseudo)
            if(j==4):
                p2=random.choice(stap)
                j=0
                    
            while (p1+"-"+p2 in pos):
                if(j==1):
                    p2=random.choice(ecoli)
                if(j==2):
                    p2=random.choice(myco)
                if(j==3):
                    p2=random.choice(pseudo)
                if(j==4):
                    p2=random.choice(stap)
                    j=0
            
            j+=1
            g.write("%s\t%s\t0\n" %(p1, p2) )
        g.close()
       
    def ppi_prediction_training(self):
        g=open("host_pathogen_interaction/training_dataset/training_dataset.txt", "w")
        for d in ['false','positive']:
            f=open("host_pathogen_interaction/training_dataset/training_"+d+".txt", "r")
            for line in f:
                g.write(line)
            f.close()
        
        params={}
        params['name']="Experiment for context-based PPI network construction"
        params['description']="experiment"
        params['email']="ycfrenchgirl2@gmail.com"
        params['owner']="Experiment"
        ds={}
        ds["prefix"]="training_"
        ds["folder"]="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/host_pathogen_interaction/training_dataset/"
        ds["description"]="training host pathogen dataset for models, labels were randomly selected"
        params["datasets"]=[ds]
        with open('predprin/params_hostpathogen.json', 'w') as fp:
            json.dump(params, fp)
            
        f=open("execute_predppi.sh","w")
        f.write("#!/bin/bash\n")
        f.write("cd predprin/\n")
        f.write("/usr/bin/python3.8 -m luigi --module main RunPPIExperiment --parameters-file params_hostpathogen.json --mode train --model none --workers 1")
        f.close()

        os.system("rm predprin/run_experiment.txt")
        print("Executing ppi prediction with host pathogen pairs")
        os.system('bash execute_predppi.sh')
   
    def run_validation_process(self):
        f=open("run_ppivalproc.sh","w")
        f.write("cd ppipubminer/\n")
        f.write("python3 pubmed_pmc_literature_pipeline.py -em 1 -rtm1 3 -fo /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/hostpathogen_papers_valproc/ -fp training_positive_hpidb.txt -fe literature_evaluation_pairs.tsv ")
        f.close()
        os.system('bash run_ppivalproc.sh')
    
    def download_protein_info(self, folder, p):
        try:
            link = "https://www.uniprot.org/uniprot/"+p+".rdf"
            f = urllib.request.urlopen(link)
            file = f.read()
            f=open( os.path.join(folder, "rdf_data", p+".rdf"), "w")
            f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace('>"','>'))
            f.close()
            
            new_link=""
            id_=p
            f=open(os.path.join(folder, "rdf_data", p+".rdf"),"r")
            for line in f:
                l=line.replace("\n","")
                if(l.find("replacedBy")!=-1):
                    new_link=l.split("=")[1].replace('"',"").replace("/>","")
                    break
            f.close()
            
            if(new_link!=""):
                id_=new_link.split("/")[-1]
                f = urllib.request.urlopen(new_link+".rdf")
                file = f.read()
                f=open(os.path.join(folder, "rdf_data", p+".rdf"),"w")
                f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace('>"','>'))
                f.close()
        except:
            print(p)
            
    def get_preferred_label(self, idu):
        label=""
        if(not os.path.isfile("predprin/rdf_data/"+idu+".rdf")):
            self.download_protein_info('predprin', idu)
            
        f=open("predprin/rdf_data/"+idu+".rdf","r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith("<skos:prefLabel")): 
                label=l.split(">")[1].replace('</skos:prefLabel','').upper()
                #print(label.upper())
        f.close()
        
        return label
        
    def build_testsDS_isolate_each_bacteria_gold(self):
        if(not os.path.isdir('host_pathogen_interaction/hpidb/')):
            #os.system("rm -rf  host_pathogen_interaction/hpidb/"+k)
            os.system("mkdir  host_pathogen_interaction/hpidb/")
            
        ppis={}
        seeds={"human": []}
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for f in folders:
            seeds[f]=[]
            
            ppis[f]=[]
            i=0
            fi=open(f+"/ppis.txt", "r")
            for line in fi:
                l=line.replace("\n","").split("\t")
                if(i>0 and (l[9]=="9606" or l[10]=='9606')): # separate only ppis whose host is human
                    id1=""
                    id2=""
                    p1=l[0].split("|")
                    for id in p1:
                        if(id.startswith("uniprotkb")):
                            id1=id.split(":")[1].split("-")[0]
                            
                    if(id1=="" and len(l)>17):
                        if(l[15].startswith("UNIPROT_AC")):
                            id1=l[15].split(":")[1]
                            
                    if(id1!=""):
                        name1=self.get_preferred_label(id1)
                        #name=l[24].split("_")[0]
                        if(name1!=''):
                            tax1=l[9]
                            if(l[9]=='9606'):
                                if(not [id1, l[9], name1] in seeds['human']):
                                    seeds['human'].append([id1, l[9], name1])
                            else:
                                if(not [id1, l[9], name1] in seeds[f]):
                                    seeds[f].append([id1, l[9], name1])  
                                          
                    p2=l[1].split("|")
                    for id in p2:
                        if(id.startswith("uniprotkb")):
                            id2=id.split(":")[1].split("-")[0]
                            
                    if(id2=="" and len(l)>17):
                        if(l[16].startswith("UNIPROT_AC")):
                            id2=l[16].split(":")[1]
                            
                    if(id2!=""):
                        name2=self.get_preferred_label(id2)
                        #name=l[25].split("_")[0]
                        if(name2!=''):
                            tax2=l[10]
                            if(l[10]=='9606'):
                                if(not [id2, l[10], name2 ] in seeds['human']):
                                    seeds['human'].append([id2, l[10], name2])
                            else:
                                if(not [id2, l[10], name2 ] in seeds[f]):
                                    seeds[f].append([id2, l[10], name2])
                            
                    if(id1!="" and id2!="" and name1!="" and name2!="" and (tax1=="9606" or tax2=='9606')):
                        ppis[f].append([id1, id2, tax1, tax2])
                        
                i+=1
            fi.close()
            
        # Prepare test datasets     
        for k in ppis.keys():
            if(not os.path.isdir('host_pathogen_interaction/hpidb/'+k)):
                #os.system("rm -rf  host_pathogen_interaction/hpidb/"+k)
                os.system("mkdir  host_pathogen_interaction/hpidb/"+k)
            
            
            print('gold---', k, len(ppis[k]))
            cut = round(len(ppis[k])/2)
            
            fx=open("host_pathogen_interaction/hpidb/"+k+"/complete_dataset_ppi.tsv","w")
                
            f=open("host_pathogen_interaction/hpidb/"+k+"/hpidb_dataset.txt","w")
            for i in range(cut):
                fx.write("1\n")
                
                with open("host_pathogen_interaction/hpidb/"+k+"/hpidb_positive.txt", "a") as fg:
                    fg.write("%s\t%s\t%s\t%s\t1\n" %(ppis[k][i][0], ppis[k][i][1], ppis[k][i][2], ppis[k][i][3]) )
                    f.write("%s\t%s\t%s\t%s\t1\n" %(ppis[k][i][0], ppis[k][i][1], ppis[k][i][2], ppis[k][i][3]) )
                    
            for i in range(cut, len(ppis[k]) ):
                fx.write("1\n")
                
                with open("host_pathogen_interaction/hpidb/"+k+"/hpidb_false.txt", "a") as fg:
                    fg.write("%s\t%s\t%s\t%s\t0\n" %(ppis[k][i][0], ppis[k][i][1], ppis[k][i][2], ppis[k][i][3]) )
                    f.write("%s\t%s\t%s\t%s\t0\n" %(ppis[k][i][0], ppis[k][i][1], ppis[k][i][2], ppis[k][i][3]) )
            f.close()
            
            fx.close()
        
        params={}
        params['name']="Experiment for context-based PPI network construction"
        params['description']="experiment string bacteria "
        params['email']="ycfrenchgirl2@gmail.com"
        params['owner']="Yasmmin"   
        datasets=[]
        # Prepare running configuration
        for k in ppis.keys():
            ds={}
            ds['name'] = "hpidb "+k
            ds["prefix"]="hpidb_"
            ds["folder"]="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/host_pathogen_interaction/hpidb/"+k+"/"
            ds["description"]="training host pathogen dataset for models, labels were randomly selected"
            datasets.append(ds)
            
        params["datasets"]=datasets
        with open('host_pathogen_interaction/hpidb/hp_params_hpidb_triplificator.json', 'w') as fp:
            json.dump(params, fp)
        os.system('cp host_pathogen_interaction/hpidb/hp_params_hpidb_triplificator.json ppi_triplification_process/hp_params_hpidb_triplificator.json')
        
        # Prepare seeds for network enrichment string
        if(not os.path.isdir('seeds')):
            os.system("mkdir seeds")
        for k in seeds.keys():
            taxon=seeds[k][0][1] 
            i=0
            f=open("seeds/"+k+"-"+taxon+".tsv", "w")
            for s in seeds[k]:
                #print(s)
                f.write("%s\t%s\n" %( s[0], seeds[k][i][2] ) )
                i+=1
            f.close()
            
    def get_string_network_from_seeds(self):
        params={}
        params['name']="Experiment for context-based PPI network construction string "
        params['description']="experiment string enriched from species seeds"
        params['email']="ycfrenchgirl2@gmail.com"
        params['owner']="Yasmmin"
        datasets=[]
        
        if(not os.path.isdir('string')):
            os.system('mkdir string')
            
        folders=['human', 'ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for f in folders:
            if(not os.path.isdir('string/'+f)):
                os.system('mkdir string/'+f)
            
            ds={}
            ds['name'] = "string_"+f
            ds["prefix"]="string_"
            ds["folder"]="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/string/"+f+"/"
            ds["description"]="test individual ppis "
            datasets.append(ds)
            
            mapp={}
            seeds=[]
            for ss in os.listdir("seeds"):
                if(ss.startswith(f)):
                    k=ss
                    taxid=k.split(".")[0].split("-")[1]
                    
            fi=open("seeds/"+k, "r")
            for line in fi:
                l=line.replace("\n","").split("\t")
                seeds.append(l[1])
            fi.close()
            
            fi=open(f+"/mapp.txt", "r")
            for line in fi:
                l=line.replace("\n","").split("\t")
                if(l[2].lower().find("uniprot_ac")!=-1):
                    mapp[l[0]]=l[1]
            fi.close()
            
            string_api_url = "https://string-db.org/api"
            output_format = "tsv-no-header"
            method = "interaction_partners"
            request_url = "/".join([string_api_url, output_format, method])
            paramsr = {
               "identifiers" : "%0d".join(seeds), # your protein
                "species" : int(taxid), # species NCBI identifier 
                "limit" : 100,
                "required_score" : 900,
                "caller_identity" : "www.exptodrugrep.org" # your app name
            }
            response = requests.post(request_url, data=paramsr)
            
            os.system('cp '+f+'/literature_link.tsv string/'+f+'/literature_link.tsv')
            
            x=open('string/'+f+"/complete_dataset_ppi.tsv","w")
            y=open('string/'+f+"/string_dataset.txt","w")
            z=open('string/'+f+"/tm_input.txt","w")
            
            n=0
            for line in response.text.strip().split("\n"):
                l = line.strip().split("\t")
                strv=l[5:]+[l[4]]
                if(l[0] in mapp.keys() and l[1] in mapp.keys()):
                    x.write('\t'.join(strv)+"\n")
                    y.write("%s\t%s\t0\n" %(mapp[l[0]], mapp[l[1]]) )
                    z.write("%s\t%s\n" %(l[2], l[3]) )
                    n+=1
            print("------------",f+":", n)
            x.close()
            y.close()
            z.close()
            
        params["datasets"]=datasets
        with open('string/hp_params_string_triplificator.json', 'w') as fp:
            json.dump(params, fp)
        os.system('cp string/hp_params_string_triplificator.json ppi_triplification_process/hp_params_string_triplificator.json')
            
    def run_ppipubminer_hostpathogen_gold(self):
        folders=['human', 'ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for f in folders:
            print(f)
            os.system('bash run_ppivalproc.sh 0 '+f+'/ tm_input.txt')
    
    def get_other_ids(self, fo, pat, mapp):
        keys=[]
        f=open(fo+"/string_dataset.txt","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            keys.append(l[0])
            keys.append(l[1])
        f.close()
        
        values=[]
        f=open(fo+"/tm_input.txt","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            values.append(l[0])
            values.append(l[1])
            
        f.close()
        
        mapp2={}
        j=0
        for i in values:
            if(i.upper() in mapp.keys()):
                if(not keys[j] in pat):
                    pat.append(keys[j])
            j+=1
        return pat
    
    def build_new_ppis(self):
        if(os.path.isdir('new_candidates_predprin')):
            os.system('rm -rf new_candidates_predprin')
        os.system('mkdir new_candidates_predprin')
        
        params={}
        params['name']="Experiment for context-based PPI network construction hostpathogen predppi "
        params['description']="experiment predprin "
        params['email']="ycfrenchgirl2@gmail.com"
        params['owner']="Yasmmin"
        datasets=[]
        
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for fo in folders:
            pat=[]
            hum=[]
            ppis=[]
            for d in ['positive','false']:
                f=open("host_pathogen_interaction/"+fo+"_ds/training_"+d+".txt","r")
                for line in f:
                    l=line.replace("\n","").split('\t')
                    if(l[2]=='9606'):
                        if(not l[0] in hum):
                            hum.append(l[0])
                        if(not l[1] in pat):
                            pat.append(l[1])
                        ppis.append([l[0], l[1]])
                    else:
                        if(not l[1] in hum):
                            hum.append(l[1])
                        if(not l[0] in pat):
                            pat.append(l[0])
                        ppis.append([l[1], l[0]])
                f.close()
            
            for ss in os.listdir("seeds"):
                if(ss.startswith(fo)):
                    k=ss
            mapp={}        
            fi=open("seeds/"+k, "r")
            for line in fi:
                l=line.replace("\n","").split("\t")
                mapp[l[1]]=l[0]
            fi.close()
            
            pat=self.get_other_ids(fo, pat, mapp)
            #print(pat)
            interactors_pat={}
            f=open(fo+"/string_dataset.txt","r")
            for line in f:
                l=line.replace('\n','').split('\t')
                if(l[0] in pat):
                    if(not l[0] in interactors_pat.keys()):
                        interactors_pat[l[0]]=[]
                    if(not l[1] in interactors_pat[l[0]]):
                        interactors_pat[l[0]].append(l[1])
                        
                if(l[1] in pat):
                    if(not l[1] in interactors_pat.keys()):
                        interactors_pat[l[1]]=[]
                    if(not l[0] in interactors_pat[l[1]]):
                        interactors_pat[l[1]].append(l[0])
            f.close()
            
            mapp={}        
            fi=open("seeds/human-9606.tsv", "r")
            for line in fi:
                l=line.replace("\n","").split("\t")
                mapp[l[1]]=l[0]
            fi.close()
            hum=self.get_other_ids('human', hum, mapp)
            
            interactors_hum={}
            f=open("human/string_dataset.txt","r")
            for line in f:
                l=line.replace('\n','').split('\t')
                if(l[0] in hum):
                    if(not l[0] in interactors_hum.keys()):
                        interactors_hum[l[0]]=[]
                    if(not l[1] in interactors_hum[l[0]]):
                        interactors_hum[l[0]].append(l[1])
                        
                if(l[1] in hum):
                    if(not l[1] in interactors_hum.keys()):
                        interactors_hum[l[1]]=[]
                    if(not l[0] in interactors_hum[l[1]]):
                        interactors_hum[l[1]].append(l[0])
            f.close()
            
            print('-------', fo)
            #for k in interactors_pat.keys():
            #    print('\tpathogen----'+k, len(interactors_pat))
            #for w in interactors_hum.keys():
            #    print('\thuman----'+w, len(interactors_hum))
            if(not os.path.isdir('new_candidates_predprin/'+fo)):
                os.system('mkdir new_candidates_predprin/'+fo)
                    
            with open('new_candidates_predprin/'+fo+"/pathogen_partners.json", 'w') as fp:
                json.dump(interactors_pat, fp)
            with open('new_candidates_predprin/'+fo+"/host_partners.json", 'w') as fp:
                json.dump(interactors_hum, fp)
                
            #print(pat)
            #print(hum)
            #print(interactors_pat)
            #print(interactors_hum)
            
            newcand=set()
            for k in interactors_pat.keys():
                for w in interactors_hum.keys():
                    #if(not [k,w] in ppis and not k+"_"+w in newcand):
                    #    newcand.add(k+"_"+w)
                    for i in interactors_pat[k]:
                        for j in interactors_hum[w]:
                            if(not [i,j] in ppis and not i+"_"+j in newcand and len(newcand)<250000 ):
                                newcand.add(i+"_"+j)
            if(len(newcand)>0):
                ds={}
                ds['name'] = "predprin new ppis from "+fo
                ds["prefix"]="predprin_"
                ds["folder"]="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/new_candidates_predprin/"+fo+"/"
                ds["description"]="test individual ppis "
                datasets.append(ds)
                
                print(fo, len(newcand))
                
                listnew=list(newcand)
                cut=round(len(newcand)/2)
                g=open('new_candidates_predprin/'+fo+"/predprin_dataset.txt","w")
                f=open('new_candidates_predprin/'+fo+"/predprin_positive.txt","w")
                for i in range(cut):
                    p1=listnew[i].split("_")[0]
                    p2=listnew[i].split("_")[1]
                    f.write("%s\t%s\t1\n" %(p1, p2) )
                    g.write("%s\t%s\t1\n" %(p1, p2) )
                f.close()
                
                f=open('new_candidates_predprin/'+fo+"/predprin_false.txt","w")
                for i in range(cut, len(newcand)):
                    p1=listnew[i].split("_")[0]
                    p2=listnew[i].split("_")[1]
                    f.write("%s\t%s\t0\n" %(p1, p2) )
                    g.write("%s\t%s\t0\n" %(p1, p2) )
                f.close()
                g.close()
                
                #os.system('cp v1_new_candidates/'+fo+'/predictions.tsv new_candidates_predprin/'+fo+"/predictions.tsv")
                #os.system('cp v1_new_candidates/'+fo+'/dataset_ppi.txt new_candidates_predprin/'+fo+"/dataset_ppi.txt")
                
        params["datasets"]=datasets
        with open('new_candidates_predprin/hp_params_predprin_triplificator.json', 'w') as fp:
            json.dump(params, fp)
        os.system('cp new_candidates_predprin/hp_params_predprin_triplificator.json ppi_triplification_process/params_predprin_triplificator.json')
         
    def mapp_new_candidates_to_hgnc_ppipubminer(self):
        "Maaping to hgnc"
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        #folders=[ 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for fo in folders:
            not_=[]
            ok=[]
            total=0
            
            print('---',fo)
            mapp={}
            x=open("new_candidates_predprin/"+fo+"/tm_input.txt","w")
            f=open("new_candidates_predprin/"+fo+"/predictions.tsv","r")
            for line in f:
                l=line.replace('\n','').split('\t')
                if(float(l[2])>0.9):
                    total+=1
                    n1=self.get_preferred_label(l[0].split('-')[0])
                    n2=self.get_preferred_label(l[1].split('-')[0])
                    if(n1!="" and n2!=""):
                        mapp[l[0].split('-')[0]+"|"+l[1].split('-')[0]]=n1+"|"+n2
                        x.write('%s\t%s\t1\n' %(n1, n2) )
                        ok.append( [l[0], l[1]] )
                    else:
                        not_.append( [l[0], l[1]] )
            f.close()
            x.close()
            
            print('not ok', len(not_))
            print('ok', len(ok))
            print('total', total)
        
            f=open("new_candidates_predprin/"+fo+"/mapping.txt", "w")
            for k in mapp.keys():
                f.write('%s\t%s\n' %(k, mapp[k]))
            f.close()
   
    def run_ppipubminer_new_candidates(self):
        "Running ppipubminer"
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        #folders=['ecoli', 'mycobacterium_tuberculosis']
        for fo in folders:
            print('---',fo)
            os.system('bash run_ppivalproc.sh 0 new_candidates_predprin/'+fo+'/ tm_input.txt')
    
    def build_literature_link_candidates(self):
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        #folders=['ecoli', 'mycobacterium_tuberculosis']
        for fo in folders:
            print('---',fo)
            mapp={}
            f=open('new_candidates_predprin/'+fo+"/mapping.txt",'r')
            for line in f:
                l=line.replace('\n','').split('\t')
                mapp[l[1].replace('|','-').lower()]=l[0].split('|')
            f.close()
            
            ok_pairs=0
            g=open('new_candidates_predprin/'+fo+"/literature_link.tsv","w")
            for f in os.listdir('new_candidates_predprin/'+fo+'/processed_sentences'):
                pair=f.split('.')[0].replace('scs_','')
                if(pair in mapp.keys()):
                    articles=[]
                    c=0
                    fg=open('new_candidates_predprin/'+fo+'/processed_sentences/'+f,'r')
                    for line2 in fg:
                        l2=line2.replace('\n','').split('\t')
                        if(not l2[0] in articles and c>0):
                            articles.append(l2[0])
                        c+=1
                    fg.close()
                    
                    uniprot=mapp[pair]
                    g.write('%s\t%s\t%s\n' %(uniprot[0], uniprot[1], ' '.join(articles)) )
                ok_pairs+=1
            g.close()
            
            if(ok_pairs==0):
                os.system('rm new_candidates_predprin/'+fo+"/literature_link.tsv")
         
    def mapp_proteins_to_hgnc_ppipubminer(self):
        mapp={}
        x=open("gold_all+hostpathogen_papers_valproc/tm_input.txt","w")
        f=open("gold_all+hostpathogen_papers_valproc/training_positive_hpidb.txt","r")
        for line in f:
            l=line.replace('\n','').split('\t')
            n1=self.get_preferred_label(l[0].split('-')[0])
            n2=self.get_preferred_label(l[1].split('-')[0])
            if(n1!="" and n2!=""):
                mapp[l[0].split('-')[0]]=n1
                mapp[l[1].split('-')[0]]=n2
                x.write('%s\t%s\t1\n' %(n1, n2) )
            else:
                print(l[0], l[1])
        f.close()
        x.close()
        
        f=open("gold_all+hostpathogen_papers_valproc/mapping.txt", "w")
        for k in mapp.keys():
            f.write('%s\t%s\n' %(k, mapp[k]))
        f.close()
   
    def run_ppipubminer_gold(self):
        os.system('bash run_ppivalproc.sh 0 gold_all+hostpathogen_papers_valproc/ tm_input.txt')
        
    def check_string_tm_score(self):
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        #folders=['ecoli', 'mycobacterium_tuberculosis']
        for fo in folders:
            print('---',fo)
            g=open(fo+"/literature_link.tsv","w")
            h=open(fo+"/tmscore_matching_study.tsv","w")
            
            tmvalues=[]
            f=open(fo+"/complete_dataset_ppi.tsv","r")
            for line in f:
                l=line.replace("\n","").split('\t')
                tmvalues.append( float(l[6]) )
            f.close()
            
            keys=[]
            f=open(fo+"/string_dataset.txt","r")
            for line in f:
                l=line.replace("\n","").split('\t')
                keys.append(l[0]+"-"+l[1])
            f.close()
            
            values=[]
            f=open(fo+"/tm_input.txt","r")
            for line in f:
                l=line.replace("\n","").split('\t')
                values.append(l[0].lower()+"-"+l[1].lower())
            f.close()
            
            with_={}
            for passed in os.listdir(fo+'/processed_sentences'):
                name=passed.split('.')[0].replace('scs_','')
                ind=values.index(name)
                pair=keys[ind]
                p1=pair.split('-')[0]
                p2=pair.split('-')[1]
                
                articles=[]
                i=0
                fg=open(fo+'/processed_sentences/'+passed,'r')
                for line in fg:
                    if(i>0):
                        l=line.replace('\n','').split('\t')
                        if(not l[0] in articles):
                            articles.append(l[0])
                    i+=1
                fg.close()
                g.write("%s\t%s\t%s\n" %(p1, p2, ' '.join(articles) ) )
                with_[pair]=len(articles)
                
            j=0
            res={'both': 0, 'ppipubminer_high': 0, 'tmscore_high': 0, 'none': 0 }
            for pair in keys:
                p1=pair.split('-')[0]
                p2=pair.split('-')[1]
                
                ppipubminer_score=0
                if(pair in with_.keys()):
                    if(with_[pair]>0):
                        ppipubminer_score=1
                        
                if(tmvalues[j]>=0.8):
                    if(ppipubminer_score==1):
                        res['both']+=1
                    else:
                        res['tmscore_high']+=1
                else:
                    if(ppipubminer_score==1):
                        res['ppipubminer_high']+=1
                    else:
                        res['none']+=1
                        
                h.write("%s\t%s\t%d\t%.2f\n" %(p1, p2, ppipubminer_score, tmvalues[j]) )
                j+=1
            
            print('Summary')
            for k in res.keys():
                print('\t'+k, res[k])
            g.close()
            h.close()
            
    def run_triplification_preparation(self):
        mapp={'predprin': '1', 'string': '2', 'hint': '3', 'hpidb': '4'}
        for f in os.listdir('ppi_triplification_process'):
            if(f.startswith('hp_') and f.find('predprin')!=-1):
                key=f.split('_')[2]
                os.system('python3 ppi_triplification_process/prepare_data_triplification.py -rt '+mapp[key]+' -fec ppi_triplification_process/'+f+' -ant /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/predprin/annotation_data/ ')
          
    def download_sequence_info(self, folder, p):
        try:
            link = "https://www.uniprot.org/uniprot/"+p+".fasta"
            f = urllib.request.urlopen(link)
            file = f.read()
            f=open("tempseq.txt","w")
            f.writelines(str(file))
            f.close()
            
            f=open( os.path.join(folder,"sequence_data", p+".fasta"), "w")
            f.close()

            f=open("tempseq.txt","r")
            for line in f:
                with open( os.path.join(folder,"sequence_data", p+".fasta"), "a") as myfile:
                    myfile.write(">"+p+"\n")
                l=str(line).replace("b'","").replace("'","").split("\\n")
                c=0
                for l_ in l:
                    if(l_!="" and l_.find(">")==-1):
                        with open( os.path.join(folder,"sequence_data", p+".fasta"), "a") as myfile:
                            myfile.write(l_)
                    c+=1
                    if(c==len(l)):
                        with open( os.path.join(folder,"sequence_data", p+".fasta"), "a") as myfile:
                            myfile.write("\n")
            f.close()
        except:
            print(p)
            pass
            
    def get_missing_sequences(self):
        proteins=set()
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        #folders=['ecoli', 'mycobacterium_tuberculosis']
        for fo in folders:
            f=open('new_candidates_predprin/'+fo+'/predprin_dataset.txt','r')
            for line in f:
                l=line.split('\t')
                if(not l[0] in proteins and not os.path.isfile('predprin/sequence_data/'+l[0]+'.fasta')):
                    proteins.add(l[0])
                if(not l[1] in proteins and not os.path.isfile('predprin/sequence_data/'+l[1]+'.fasta')):
                    proteins.add(l[1])
                    
            f.close()
        c=1
        for p in proteins:
            print(c, '/', len(proteins))
            self.download_sequence_info('predprin/', p)
            c+=1
    
    def make_info_file_proteins(self):
        folders=['ecoli', 'mycobacterium_tuberculosis', 'staphylococcus_aureus','pseudomonas_aeruginosa']
        for fo in folders:
            proteins=set()
            f=open("host_pathogen_interaction/hpidb/"+fo+"/hpidb_dataset.txt", "r")
            for line in f:
                l=line.replace('\n','').split('\t')
                if(not l[0] in proteins):
                    proteins.add(l[0])
                if(not l[1] in proteins):
                    proteins.add(l[1])
            f.close()
            
            annotation_path='predprin/annotation_data/'
            f=open("host_pathogen_interaction/hpidb/"+fo+"/info_proteins.tsv","w")
            for pr in proteins:
                if(os.path.isfile(annotation_path+pr+".tsv")):
                    g=open(annotation_path+pr+".tsv","r")
                    for line in g:
                        l=line.replace("\n","")
                        if(line!=""):
                            f.write(line)
                    g.close()
            f.close()
         
# python3 -m luigi --module main RunPPIExperiment --parameters-file hp_params_predprin_triplificator.json --mode test --model /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/portfolio/bioinformatics_advances/host_pathogen_interaction/training_dataset/model_trained.joblib --workers 4
            
a=HostPathogenBacteria()
#a.build_mapping_training_data()
#a.complement_positive()
#a.build_false_dataset()
#a.ppi_prediction_training()
#a.run_validation_process()

# getting models
#a.build_train_individual_datasets()
"""
-----ecoli 277
-----mycobacterium_tuberculosis 22
-----staphylococcus_aureus 10155
-----pseudomonas_aeruginosa 31
"""
#a.run_train_ppi()

#a.build_testsDS_isolate_each_bacteria_gold() # filtering human host hpidb

"""
gold--- ecoli 191
gold--- mycobacterium_tuberculosis 15
gold--- staphylococcus_aureus 39
gold--- pseudomonas_aeruginosa 18
"""

a.make_info_file_proteins()

#a.get_string_network_from_seeds()
"""
------------ human: 5142
------------ ecoli: 146
------------ mycobacterium_tuberculosis: 48
------------ staphylococcus_aureus: 4
------------ pseudomonas_aeruginosa: 1
"""
#a.check_string_tm_score()

#a.run_ppipubminer_hostpathogen_gold()

#a.run_ppipubminer_gold()

#a.mapp_proteins_to_hgnc_ppipubminer()

#a.build_new_ppis()
"""
ecoli 230547
mycobacterium_tuberculosis 66933
staphylococcus_aureus 9916
pseudomonas_aeruginosa 2479
"""


#a.get_missing_sequences()

a.run_triplification_preparation()
#a.build_literature_link_candidates()

#a.mapp_new_candidates_to_hgnc_ppipubminer()
#a.run_ppipubminer_new_candidates()
"""
--- ecoli
not ok 32
ok 6155
total 6187
--- mycobacterium_tuberculosis
not ok 2
ok 4205
total 4207
--- staphylococcus_aureus
not ok 61
ok 55
total 116
--- pseudomonas_aeruginosa
not ok 0
ok 87
total 87

"""


