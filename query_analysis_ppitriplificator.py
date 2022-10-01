from franz.openrdf.connect import ag_connect
from franz.openrdf.query.query import QueryLanguage

class Query_scenarios:
    conn=None
    def __init__(self):
        self.conn = ag_connect('ppitriplificator', create=False, clear=False)

    def run_query_1_ds_organisms(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            select ?organism (group_concat(?name; separator=" | ") as ?names) where {
                ?dataset a ontoppi:PPIDataset . 
                ?dataset ontoppi:fromOrganism ?organism . 
                ?dataset rdfs:label ?name . 
            } group by ?organism
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s" % (binding_set.getValue("organism"), binding_set.getValue("names") ))

    def run_query_2_interactions_literature(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            select ?interaction ?decision ?pubmed where {
                ?interaction ontoppi:hasScore ?decision .
                ?decision ontoppi:supportingArticle ?pubmed . 
            } 
        """
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            select ?interaction ?uniprot1 ?uniprot2 (group_concat(?pubmed; separator=" | ") as ?pubmeds) where {
                ?interaction ontoppi:hasScore ?decision . 
                ?decision ontoppi:supportingArticle ?pubmed . 
                
                ?interaction  ontoppi:participant1 ?protein1 . 
                ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 

                ?interaction  ontoppi:participant2 ?protein2 . 
                ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
            
            } group by ?interaction ?uniprot1 ?uniprot2
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s %s %s" % (binding_set.getValue("interaction"), binding_set.getValue("uniprot1"), binding_set.getValue("uniprot2"), binding_set.getValue("pubmeds") ) )

    def run_query_3_0_order_bps_bacteria(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select distinct ?interaction ?protein where {
                ?dataset ontoppi:fromOrganism ppiprov:organism_escherichia_coli . 
                ?interaction ontoppi:belongsTo ?dataset . 
                ?interaction ontoppi:participant1 ?protein
            }  limit 5
        """
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select distinct ?bp (count(?bp) as ?count) where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_escherichia_coli  . 
                ?interaction  ontoppi:participant1 ?protein . 
                ?protein ontoppi:hasGO_bp_annotation ?bp . 
            } group by ?bp order by desc(?count) limit 5
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s" % (binding_set.getValue("bp"), binding_set.getValue("count") ) )

    def run_query_3_interactions_specific_bp_bacteria(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select ?interaction ?uniprot1 ?uniprot2 where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_escherichia_coli  . 

                ?interaction  ontoppi:participant1 ?protein1 . 
                ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                ?protein1 ontoppi:hasGO_bp_annotation ?bp1 . 

                ?interaction  ontoppi:participant2 ?protein2 . 
                ?protein2 ontoppi:hasGO_bp_annotation ?bp2 . 
                ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
            
                filter (regex(?bp1, 'GO:0006355', 'i' ) &&  regex(?bp2, 'GO:0006355', 'i' )) .
            } 
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s %s" % (binding_set.getValue("interaction"), binding_set.getValue("uniprot1"), binding_set.getValue("uniprot2") ) )

    def run_query_4_interactions_scores_specificBP_bacteria(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

            select distinct ?uniprot1 ?uniprot2 (group_concat( ?method_name; separator=" | ") as ?methods) (group_concat( ?value; separator=" | ") as ?scores) where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_escherichia_coli . 

                ?interaction ontoppi:hasScore ?decision . 
                ?decision ontoppi:predictedValue ?value . 
                ?decision ontoppi:basedOn ?method . 
                ?method rdfs:label ?method_name . 
                
                ?interaction  ontoppi:participant1 ?protein1 . 
                ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                ?protein1 ontoppi:hasGO_bp_annotation ?bp1 . 

                ?interaction  ontoppi:participant2 ?protein2 . 
                ?protein2 ontoppi:hasGO_bp_annotation ?bp2 . 
                ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
            
                filter (regex(?bp1, 'GO:0006355', 'i' ) &&  regex(?bp2, 'GO:0006355', 'i' )) .
            } group by ?uniprot1 ?uniprot2
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        tuple_query.setIncludeInferred(True) #if false it does not retrieve inferred triples
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s %s %s" % (binding_set.getValue("uniprot1"), binding_set.getValue("uniprot2"), binding_set.getValue("methods"), binding_set.getValue("scores") ) )  

    def run_query_5_0_order_bps_human(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select distinct ?interaction ?protein where {
                ?dataset ontoppi:fromOrganism ppiprov:organism_escherichia_coli . 
                ?interaction ontoppi:belongsTo ?dataset . 
                ?interaction ontoppi:participant1 ?protein
            }  limit 5
        """
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select distinct ?bp (count(?bp) as ?count) where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_homo_sapiens  . 
                ?interaction  ontoppi:participant1 ?protein . 
                ?protein ontoppi:hasGO_bp_annotation ?bp . 
            } group by ?bp order by desc(?count) limit 5
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s" % (binding_set.getValue("bp"), binding_set.getValue("count") ) )

    def run_query_5_interactions_specific_bp_human(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>

            select ?interaction ?uniprot1 ?uniprot2 where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_homo_sapiens  . 

                ?interaction  ontoppi:participant1 ?protein1 . 
                ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                ?protein1 ontoppi:hasGO_bp_annotation ?bp1 . 

                ?interaction  ontoppi:participant2 ?protein2 . 
                ?protein2 ontoppi:hasGO_bp_annotation ?bp2 . 
                ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
            
                filter (regex(?bp1, 'GO:0045944', 'i' ) &&  regex(?bp2, 'GO:0045944', 'i' )) .
            } 
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s %s" % (binding_set.getValue("interaction"), binding_set.getValue("uniprot1"), binding_set.getValue("uniprot2") ) )

    def run_query_6_interactions_scores_specificBP_human(self):
        query="""
            prefix ontoppi: <https://www.ypublish.info/protein_interaction_domain_ontology#>
            prefix ppiprov: <https://www.ypublish.info/provenance_information#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>

            select distinct ?uniprot1 ?uniprot2 (group_concat(?method_name; separator=" | ") as ?methods) (group_concat(?value; separator=" | ") as ?scores) where {
                ?interaction ontoppi:belongsTo ?dataset . 
                ?dataset ontoppi:fromOrganism ppiprov:organism_homo_sapiens  . 

                ?interaction ontoppi:hasScore ?decision . 
                ?decision ontoppi:predictedValue ?value . 
                ?decision ontoppi:basedOn ?method . 
                ?method rdfs:label ?method_name . 
                
                ?interaction  ontoppi:participant1 ?protein1 . 
                ?protein1 ontoppi:hasUniprotCorrespondent ?uniprot1 . 
                ?protein1 ontoppi:hasGO_bp_annotation ?bp1 . 

                ?interaction  ontoppi:participant2 ?protein2 . 
                ?protein2 ontoppi:hasGO_bp_annotation ?bp2 . 
                ?protein2 ontoppi:hasUniprotCorrespondent ?uniprot2 . 
            
                filter (regex(?bp1, 'GO:0045944', 'i' ) &&  regex(?bp2, 'GO:0045944', 'i' )) .
            } group by ?uniprot1 ?uniprot2
        """
        tuple_query = self.conn.prepareTupleQuery(QueryLanguage.SPARQL, query)
        tuple_query.setIncludeInferred(True) #if false it does not retrieve inferred triples
        result = tuple_query.evaluate()
        with result:
            for binding_set in result:
                print("%s %s %s %s" % (binding_set.getValue("uniprot1"), binding_set.getValue("uniprot2"), binding_set.getValue("methods"), binding_set.getValue("scores") ) )        

class Running_config:
    def run(self, args):
        if(args.query_option!=""):
            if(args.query_option in [1,2,3,4,5,6,7,8]):
                q=Query_scenarios()
                
                n=args.query_option
                if(n==1):
                    q.run_query_1_ds_organisms()
                if(n==2):
                    q.run_query_2_interactions_literature()

                if(n==3):
                    q.run_query_3_0_order_bps_bacteria()
                if(n==4):
                    q.run_query_3_interactions_specific_bp_bacteria()
                if(n==5):
                    q.run_query_4_interactions_scores_specificBP_bacteria()

                if(n==6):
                    q.run_query_5_0_order_bps_human()
                if(n==7):
                    q.run_query_5_interactions_specific_bp_human()
                if(n==8):
                    q.run_query_6_interactions_scores_specificBP_human()
            else:
                print("Please, choose a valid option between 1 and 8.")
        else:
            print("Please, choose the query option you want")

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description='Tool to query PPI semantic data', formatter_class=RawTextHelpFormatter)
parser.add_argument("-q", "--query_option", action="store", help="\
1 - Get all the different organisms whose interactions are stored in the database\n\
2 - Get the interactions that have scientific papers associated and the list of these papers\n\
3 - Get a list of the most frequent biological processes annotated for the interactions of Escherichia coli bacteria\n\
4 - Get only the interactions belonging to a specific biological process (regulation of transcription, DNA-templated) in Escherichia coli bacteria\n\
5 - Get the scores of interactions belonging to a specific biological process (regulation of transcription, DNA-templated) in Escherichia coli bacteria\n\
6 - Get a list of the most frequent biological processes annotated for the interactions of human organism\n\
7 - Get only the interactions belonging to a specific biological process (positive regulation of transcription by RNA polymerase II) in human organism\n\
8 - Get the scores of interactions belonging to a specific biological process (positive regulation of transcription by RNA polymerase II) in human organism\n\
", type=int)
args = parser.parse_args()
r=Running_config()
r.run(args)


