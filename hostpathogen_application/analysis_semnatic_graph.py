from franz.openrdf.connect import ag_connect
from franz.openrdf.query.query import QueryLanguage

from plotly.subplots import make_subplots
import plotly.graph_objects as go

class ExplorationSemanticGraph:
    conn=None
    def __init__(self):
        self.conn = ag_connect('ppitriplificator', create=False, clear=False)
          
    def load_obo(self):
        goa={}
        flag=False
        f=open("go.obo","r")
        for line in f:
            l=line.replace('\n','')
            if(l.find('[Term]')!=-1):
                id_=""
                name=""
                branch=""
                flag=True
            
            if(flag):
                if(l.startswith('id: GO')):
                    id_=l.split(": ")[1]
                
                if(l.startswith('name: ')):
                    name=l.split(": ")[1]
                  
                if(l.startswith('namespace:')):
                    branch=l.split(": ")[1]
                    if(id_!="" and name!="" and branch!=""):
                        goa[id_] = { 'name': name, 'branch': branch }
                
        f.close()   
        
        return goa   
           
    def run_query_bp_hostPathogen(self):
        for ds in ['predprin', 'hpidb']:
            print('----> Dataset ', ds)
            
            mapp={ '83331': 'Mycobacterium tuberculosis', '83332': 'Mycobacterium tuberculosis', '9606': 'Homo sapiens', '83332': 'Escherichia coli', '233413':'Mycobacterium tuberculosis', '83334': 'Escherichia coli', '93061': 'Staphylococcus aureus', '233413': 'Mycobacterium tuberculosis', '208964': 'Pseudomonas aeruginosa', '287': 'Pseudomonas aeruginosa', '562': 'Escherichia coli', '1280': 'Staphylococcus aureus', '158879': 'Staphylococcus aureus', '93061': 'Staphylococcus aureus', '574521': 'Escherichia coli' }
            goa=self.load_obo()
            
            for bra in ['cc', 'bp', 'mf']:
                print('\t----> branch ', bra)
                
                res={}
                taxons=set()
                query="""
                    prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
                    prefix ppiprov: <https://www.ypublish.info/provenance_information#>
                    prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

                    select distinct ?nameds ?uniprot1 ?uniprot2 ?taxon1 ?taxon2 (group_concat(distinct ?bp1; separator=" | ") as ?processes_protein_1) (group_concat(distinct ?bp2; separator=" | ") as ?processes_protein_2) where {
                    ?interaction  ontoppi:participant1 ?protein1 . 
                    ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                    ?protein1 ontoppi:hasGO_"""+bra+"""_annotation ?bp1 . 
                    ?protein1 ontoppi:fromOrganism ?taxon1 .

                    ?interaction  ontoppi:participant2 ?protein2 . 
                    ?protein2 ontoppi:hasGO_"""+bra+"""_annotation ?bp2 . 
                    ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
                    ?protein2 ontoppi:fromOrganism ?taxon2 .
                  
                  ?interaction ontoppi:belongsTo ?dataset . 
                  ?dataset rdfs:label ?nameds .
                	
                  filter (?taxon1 != ?taxon2 && regex(?nameds, '"""+ds+"""', 'i')) .
                    
                } group by ?nameds ?uniprot1 ?uniprot2 ?taxon1 ?taxon2
                """
                tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
                tuple_query.setIncludeInferred(True) #if false it does not retrieve inferred triples
                result = tuple_query.evaluate()
                with result:
                    for binding_set in result:
                        #print("%s %s %s %s" % (str(binding_set.getValue("taxon1")), binding_set.getValue("taxon2"), binding_set.getValue("processes_protein_1"), binding_set.getValue("processes_protein_2") ) )   
                        t1 = str(binding_set.getValue("taxon1")).split('/')[-1].replace('>','')
                        t2 = str(binding_set.getValue("taxon2")).split('/')[-1].replace('>','')
                        
                        if(not t1 in taxons):
                            taxons.add(t1)
                        if(not t2 in taxons):
                            taxons.add(t2)
                        
                        if(t1=='9606' or t2=='9606'):  
                            """
                            if(t1=='9606'):
                                t1=t1+'_'+t2
                                
                            if(t2=='9606'):
                                t2=t2+'_'+t1
                                """
                            if(t1 in mapp.keys() and t2 in mapp.keys()):
                                t1=mapp[t1]
                                t2=mapp[t2]    
                                
                                if(not t1 in res.keys()):
                                    res[t1]={}
                                prs1 = str(binding_set.getValue("processes_protein_1")).replace('"','').split(' | ')
                                for bp in prs1:
                                    if(not bp in res[t1].keys() and bp in goa.keys()):
                                        res[t1][bp]=0
                                    """else:
                                        if(not bp in goa.keys()):
                                            print(t1, bp)"""
                                        
                                    if(bp in res[t1].keys()):
                                        res[t1][bp]+=1
                                        
                                if(not t2 in res.keys()):
                                    res[t2]={}
                                prs2 = str(binding_set.getValue("processes_protein_2")).replace('"','').split(' | ')
                                for bp in prs2:
                                    if(not bp in res[t2].keys() and bp in goa.keys()):
                                        res[t2][bp]=0
                                    """else:
                                        if(not bp in goa.keys()):
                                            print(t2, bp)"""
                                        
                                    if(bp in res[t2].keys()):
                                        res[t2][bp]+=1
                
                for t in taxons:
                    print(t)
                
                plts={}
                
                ins=1
                r=1
                c=1
                f=open(ds+"_"+bra+"_results_enrichment.tsv","w")
                for k in res.keys():
                    x=[]
                    y=[]       
                    sorted_=dict(sorted(res[k].items(), key=lambda item: item[1], reverse=True))
                    total=sum(res[k].values())
                    for s in sorted_:
                        perc=sorted_[s]/total
                        #print(k, s, goa[s]['name'], sorted_[s], perc)
                        if(len(x)<5):
                            x.append(goa[s]['name'])
                            y.append(perc*100)
                            
                        f.write("%s\t%s\t%s\t%i\t%.4f\n" %(k, s, goa[s]['name'], sorted_[s], perc) )
                    
                    """if(k in mapp.keys()):
                        if(not mapp[k] in mapp.keys()):
                            plts[mapp[k]]=[x, y]"""
                           
                    plts[k]=[x, y, r, c]
                    c+=1
                    if(ins%2==0):
                        r+=1
                        c=1
                    ins+=1
                f.close()
            
            """
                labels=tuple(list(plts.keys()))
                fig = make_subplots( rows=r, cols=2, subplot_titles=labels)
                for k in plts.keys():
                    fig.add_trace(go.Bar(x=plts[k][0], y=plts[k][1]), row=plts[k][2], col=plts[k][3])    
                fig.update_layout(height=700, width=700, title_text="Most popular "+bra.upper()+" annotations")
                fig.write_image('panel_'+bra+'.png')
            """
 
"""
 prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
                    prefix ppiprov: <https://www.ypublish.info/provenance_information#>
                    prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

                    select distinct ?uniprot1 ?uniprot2 ?taxon1 ?taxon2 (group_concat(distinct ?bp1; separator=" | ") as ?processes_protein_1) (group_concat(distinct ?bp2; separator=" | ") as ?processes_protein_2) where {
                        ?interaction  ontoppi:participant1 ?protein1 . 
                        ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                        ?protein1 ontoppi:hasGO_cc_annotation ?bp1 . 
                        ?protein1 ontoppi:fromOrganism ?taxon1 .

                        ?interaction  ontoppi:participant2 ?protein2 . 
                        ?protein2 ontoppi:hasGO_cc_annotation ?bp2 . 
                        ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
                        ?protein2 ontoppi:fromOrganism ?taxon2 .
                    	
                      filter (?taxon1 != ?taxon2 ) .
                        
                    } group by ?uniprot1 ?uniprot2 ?taxon1 ?taxon2
""" 
            
"""
Distinct taxons 
prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

            select distinct ?taxon ?dataset where {
                ?protein rdf:type ontoppi:PairComponent . 
                ?protein ontoppi:fromOrganism ?taxon .
            	
            } 
"""        
"""
Experiments and datasets and the number of ppis in each of them
prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

            select  ?exp ?dataset (count(?interaction) as ?number_ppis) where {
              ?exp ontoppi:hasDataset ?dataset . 
              ?interaction ontoppi:belongsTo ?dataset .
            } group by ?exp ?dataset
"""
            
a=ExplorationSemanticGraph()
a.run_query_bp_hostPathogen()

