import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
#gbfile = "\\Users\\jasmi\\Desktop\\Bioinfo\\PKS\\pks.gbk"

path = "D:\\Studium\\A_Master\\3. Semester\\trans-AT PKS\\Datensets\\FinalesDatenset"
files = list()
for file in os.scandir(path) :
    files.append(file)
#print(files)

def retrieve_features_from_genbankfile(gbfile):
    for gbrecord in SeqIO.parse(gbfile, "genbank") :
        gbfeature = gbrecord.features       #speichern aller features in gbfeature
        anno = gbrecord.annotations
        gbcontent = gbrecord
        id_file = gbrecord.id
        
    return anno , gbfeature, gbcontent, id_file
    
def get_features(gbfeature):
    
    list_subregion_start = list()
    list_subregion_end = list()
    list_cds_in_region_start = list()
    list_cds_in_region_end = list() 
    list_cds_features = list()
    list_all_domains_start = list()
    list_all_domains_end = list()
    list_domain_features = list()
    list_KS_zero_start = list()
    list_KS_zero_end = list()
    list_ACP_beta_start = list()
    list_ACP_beta_end = list()
    list_PP_start = list()
    list_PP_end = list()
    list_module_start = list()
    list_module_end = list()
    list_modules_few_domains_len = list()
    list_modules_few_domains_start = list()
    list_modules_few_domains_end = list()
    
    for feature in gbfeature :                         #durchlaufen durch alle features
        type = feature.type
    #    print(type)
        if type == "subregion" :                         #auslesen Feature Type
            try:
                label_subregion = feature.qualifiers["label"]   #auslesen product in pks_features
                str_label_subregion = str(label_subregion)                #umwandeln product typ in string
            except:
                print("no product qualifier avaliable")
                break
            if "polyketide" or "Polyketide" in str_label_subregion :               #suche nach pk in string des labels
    #            print("PKS!!", str_product_type)
    #            print(feature)
                subregion_start = int(feature.location.start)
                if subregion_start > 0 :
                    subregion_start += 1
                list_subregion_start.append(subregion_start)
                subregion_end = int(feature.location.end)
                list_subregion_end.append(subregion_end)
        
        elif type == "CDS" :                         #auslesen Feature Type
            cds_start = int(feature.location.start)
            if cds_start > 0 :
                cds_start += 1
            cds_end = int(feature.location.end)
            list_cds_features.append(feature)
            list_cds_in_region_start.append(cds_start)
            list_cds_in_region_end.append(cds_end)
            #aa_seq_qual = feature.qualifiers["translation"]
        
        elif type == "aSDomain" :                       # Locus Bsp: [44320:44977](-)
            list_domain_features.append(feature)
            asdomain_qual = feature.qualifiers["aSDomain"]
            str_asdomain_qual = str(asdomain_qual)
            start_domain = int(feature.location.start) + 1
            end_domain = int(feature.location.end)
            list_all_domains_start.append(start_domain)
            list_all_domains_end.append(end_domain)
            
            if "PKS_KS" in str_asdomain_qual :
                aa_seq_qual_KS = feature.qualifiers["translation"]
                str_aa_seq_qual_KS = str(aa_seq_qual_KS)
                split_aa_seq_KS = str_aa_seq_qual_KS.split("'")
                str_aa_seq_KS = split_aa_seq_KS[1]
                if "HGTGT" in str_aa_seq_KS :
                    music = "music is life :) <3"
                    bio = "we are biology nerds!!"
                    trans_PKS = "where the funky domains at???"
                else :
                    list_KS_zero_start.append(start_domain)
                    list_KS_zero_end.append(end_domain)
            
            elif "ACP_beta" in str_asdomain_qual :
                list_ACP_beta_start.append(start_domain)
                list_ACP_beta_end.append(end_domain)
                
            elif "PP-binding" in str_asdomain_qual or "PKS_PP" in str_asdomain_qual:
                list_PP_start.append(start_domain)
                list_PP_end.append(end_domain)

        elif type == "aSModule":
            start_module = int(feature.location.start)
            if start_module > 0 :
                start_module = start_module + 1
            end_module = int(feature.location.end)
            list_module_start.append(start_module)
            list_module_end.append(end_module)                     # Locus Bsp: [44320:44977](-)
            asmodule_qual = feature.qualifiers["domains"]
            numb_domains = len(asmodule_qual)
            for domain in asmodule_qual :
                if ("PP-binding" in domain or "ACP" in domain or "PKS_PP" in domain) and numb_domains < 4 :
#                    print("\nFor Schleifen Domain :-) :" , domain)
#                    print(start_module, ".." , end_module , str_asmodule_qual)
                    list_modules_few_domains_len.append(numb_domains)
                    list_modules_few_domains_start.append(start_module)
                    list_modules_few_domains_end.append(end_module)
    
    return list_subregion_start, list_subregion_end, list_cds_in_region_start, list_cds_in_region_end, list_cds_features, list_all_domains_start, list_all_domains_end, list_domain_features, list_KS_zero_start, list_KS_zero_end, list_modules_few_domains_len, list_modules_few_domains_start, list_modules_few_domains_end, list_ACP_beta_start, list_ACP_beta_end, list_PP_start, list_PP_end

def get_translation_cds(gbfeature):
    list_translation_cds = list()   
    numb_cds_in_subregions = len(list_cds_in_region_start)
    count_for = 1
    cds_start = None
    for feature in gbfeature :
        count_cds_translation = 0
 #        print("FOR SCHLEIFE RUNDE:" , count_for)                    #durchlaufen durch alle features #translation zu CDS in Subregion speichern
        type = feature.type
        if type == "CDS" :                         #auslesen Feature Type
            cds_start = int(feature.location.start)
            if cds_start > 0 :
                cds_start += 1
 #            print("CDS gefunden yaaiii" , count_cds_translation , cds_start)
 #            cds_end = int(feature.location.end
            while count_cds_translation < numb_cds_in_subregions :
                if cds_start == list_cds_in_region_start[count_cds_translation] :
                    aa_seq_qual = feature.qualifiers["translation"]
                    str_aa_seq_qual = str(aa_seq_qual)
                    split_aa_seq = str_aa_seq_qual.split("'")
                    str_aa_seq = split_aa_seq[1]
                    list_translation_cds.append(str_aa_seq)
                    count_cds_translation = count_cds_translation + 1
                else :
                    count_cds_translation = count_cds_translation + 1
                #print(count_cds_translation , "+1 im eeeeelse")
            
        count_for = count_for + 1
         
    return list_translation_cds

def generate_gaps_duplicate_domains(list_domain_features):
    list_duplicate_domains_start = list()
    list_duplicate_domains_end = list()
    list_assign_domain_type = list()
    count_next_domain = 1
    for domain in list_domain_features:
        function_domain = domain.qualifiers["aSDomain"]
        str_function_domain = str(function_domain)
        split_name = str_function_domain.split("'")
        name_domain = split_name[1]
        try:
            next_domain = list_domain_features[count_next_domain]
        except:
            break
        function_next_domain = next_domain.qualifiers["aSDomain"]
        str_function_next_domain = str(function_next_domain)
        split_name_next = str_function_next_domain.split("'")
        name_next_domain = split_name_next[1]
        start_domain = int(domain.location.start) + 1
        start_next_domain = int(next_domain.location.start) + 1
        end_domain = int(domain.location.end)
        end_next_domain = int(next_domain.location.end)
        id_domain = str(start_domain)
        id_next_domain = str(start_next_domain)
        assign_type_domain = id_domain + ".duplicate_" + name_domain
        assign_type_next_domain = id_next_domain + ".duplicate_" + name_next_domain

# amp binding? condensation? 
        count_next_domain +=1
        if "Condensation" in str_function_domain or "AMP-binding" in str_function_domain:
            continue
        if "Condensation" in str_function_next_domain or "AMP-binding" in str_function_next_domain:
            continue
        if "DH" in str_function_domain:
            str_function_domain == "DH"
        if "DH" in str_function_next_domain:
            str_function_next_domain = "DH"
        if str_function_domain == str_function_next_domain:
            list_duplicate_domains_start+=[start_next_domain]
            list_duplicate_domains_end+=[end_next_domain]
            list_assign_domain_type+=[assign_type_next_domain]
            #print("______DOUBLEDOMAIN______")
            #print(f"Runde {count_next_domain}: current domain: {start_domain} {str_function_domain}, next domain: {start_next_domain} {str_function_next_domain}")    

    return list_duplicate_domains_start, list_duplicate_domains_end, list_assign_domain_type

def cut_domains_longer_than_one_cds(list_cds_in_region_start, list_cds_in_region_end, list_gaps_start, list_gaps_end, list_gaps_type):
    count_gaps = 0
    list_cut_domains_end = list()
    list_cut_domains_start = list()
    for gap_start in list_gaps_start :
        count_cds = 0
        numb_cds = int(len(list_cds_in_region_start))
        while count_cds < numb_cds :
            cds_start = list_cds_in_region_start[count_cds]
            cds_end = list_cds_in_region_end[count_cds]
            gap_type = list_gaps_type[count_gaps]
            gap_end = list_gaps_end[count_gaps]
            if "before" in gap_type:
                if cds_start <= gap_end <= cds_end:
                    if gap_start < cds_start:
                        gap_start = cds_start
                    list_cut_domains_start.append(gap_start)
                    list_cut_domains_end.append(gap_end)
                    count_cds = numb_cds
                else:
                    count_cds += 1
            elif "after" in gap_type:
                if cds_start <= gap_start <= cds_end:
                    if gap_end > cds_end:
                        gap_end = cds_end
                    list_cut_domains_start.append(gap_start)
                    list_cut_domains_end.append(gap_end)
                    count_cds = numb_cds
                else:
                    count_cds += 1
            else:
                list_cut_domains_start.append(gap_start)
                list_cut_domains_end.append(gap_end)
                count_cds = numb_cds                    
                
        count_gaps += 1
        
    return list_cut_domains_start, list_cut_domains_end
                
def generate_gaps_before_and_after_domains(list_specific_domain_start, list_all_domains_start, list_all_domains_end, list_domain_features):
    list_all_gaps_start = list()          #Listen in die später alle erstellten Lücken eingelesen werden
    list_all_gaps_end = list()
    list_domains_for_gaps_start = list()
    list_domains_for_gaps_end = list()
    count_domains = 0
    count_domains_before = -1
    count_domains_after = 1
    numb_specific_domains = len(list_specific_domain_start)
    numb_domains_in_region = len(list_all_domains_start)
    list_domain_type = list()
    list_assign_domain_type = list()
    for domain_start in list_specific_domain_start:
        count_while = 0
        len_domain_list = len(list_domains_for_gaps_start)
        if len_domain_list > 2 :                        #Wenn liste 3 items (?) hat (KS zero + domäne davor & danach)
            list_domains_for_gaps_start = list()        #und somit die lücken erstellt wurden wird liste geleert
            list_domains_for_gaps_end = list()
            len_domain_list = len(list_domains_for_gaps_start)
#        print(KS_zero_start , list_all_domains_start[count_domains])
        while count_while < numb_specific_domains :
#            print( "KS zero, Start matching domain:" , KS_zero_start , list_all_domains_start[count_domains])
            current_domain_start = list_all_domains_start[count_domains]
            current_domain = list_domain_features[count_domains]
            domain_type = current_domain.qualifiers["aSDomain"]
            str_domain_type = str(domain_type)
            if "PP" in str_domain_type:
                str_domain_type = "PP-binding"
            if "PKS_KS" in str_domain_type:
                str_domain_type = "Non_Elongating_KS"
            if "ACP_beta" in str_domain_type:
                str_domain_type = "ACP_beta"
            if domain_start == current_domain_start:

                type_before = ".before_" + str_domain_type
                type_after = ".after_" + str_domain_type
#                print("In the if!!")
#                print( "KS zero, Start matching domain:" , KS_zero_start , list_all_domains_start[count_domains])
                if count_domains_after < numb_domains_in_region:
                    #s_before, s, s_after = list_all_domains_start[count_domains_before], list_all_domains_start[count_domains], list_all_domains_start[count_domains_after]
                    #e_before, e, e_after = list_all_domains_end[count_domains_before], list_all_domains_end[count_domains], list_all_domains_end[count_domains_after]
                    list_domains_for_gaps_start+=[list_all_domains_start[count_domains_before],list_all_domains_start[count_domains],list_all_domains_start[count_domains_after]]
                    list_domains_for_gaps_end+=[list_all_domains_end[count_domains_before],list_all_domains_end[count_domains],list_all_domains_end[count_domains_after]]
                    list_domain_type+=[type_before,type_after]
                    count_while = numb_specific_domains
                    len_domain_list = len(list_domains_for_gaps_start)
                elif count_domains_after >= numb_domains_in_region:
                    list_domains_for_gaps_start+=[list_all_domains_start[count_domains_before],list_all_domains_start[count_domains]]
                    list_domains_for_gaps_end+=[list_all_domains_end[count_domains_before],list_all_domains_end[count_domains]]
                    list_domain_type+=[type_before]
                    count_while = numb_specific_domains
                    len_domain_list = (len(list_domains_for_gaps_start) + 1)
                elif count_domains_before < 0 :
                    list_domains_for_gaps_start+=[list_all_domains_start[count_domains],list_all_domains_start[count_domains_after]]
                    list_domains_for_gaps_end+=[list_all_domains_end[count_domains],list_all_domains_end[count_domains_after]]
                    list_domain_type+=[type_after]
                    count_while = numb_specific_domains
                    len_domain_list = (len(list_domains_for_gaps_start) + 1)
            else :
#                print("Else, everything + ooooone")
                count_domains += 1
                count_domains_before += 1
                count_domains_after += 1

            if len_domain_list > 2 :
                count_end_gap = 1
                count_type = 0
                for domain_end in list_domains_for_gaps_end :
                    try :
                        gap_end = list_domains_for_gaps_start[count_end_gap] - 1
                    except :
                        continue
                    gap_start = domain_end
                    id_domain = str(gap_start)
                    add_type = id_domain+list_domain_type[count_type]
                    list_assign_domain_type.append(add_type)
                    list_all_gaps_start.append(gap_start)
                    list_all_gaps_end.append(gap_end)
#                    print("Ende vorherhiger Domain & Start Nächster:    " , domain_end , list_KS_zero_domains_for_gaps_start[count_end_gap])
 #                   print("Start & Ende Gap:                       " , gap_start , gap_end)
                    count_end_gap += 1
                    count_type += 1
                    
    return list_all_gaps_start, list_all_gaps_end, list_assign_domain_type

def generate_small_module_gaps(list_modules_few_domains_len, list_modules_few_domains_start, list_modules_few_domains_end, list_all_domains_start, list_all_domains_end):
    list_domain_small_module_start = list()
    list_domain_small_module_end = list()
    list_all_gaps_start = list()
    list_all_gaps_end = list()
    list_assign_domain_type = list()
    count_domain_end = 0
    count_numb_domains = 0
    count_id_module = 1
    for domain_start in list_all_domains_start :
        domains_added = len(list_domain_small_module_start)
        count_module = 0
        domain_end = list_all_domains_end[count_domain_end]
        try :
            numb_domains_in_module = list_modules_few_domains_len[count_numb_domains]
        except :
            break
        if domains_added == numb_domains_in_module:
            list_domain_small_module_start = list()
            list_domain_small_module_end = list()
            count_numb_domains += 1
            count_id_module += 1
            domains_added = len(list_domain_small_module_start)
        try :
            numb_domains_in_module = list_modules_few_domains_len[count_numb_domains]
            numb_modules = list_modules_few_domains_start[count_module]
        except :
            break
        numb_modules = len(list_modules_few_domains_start)
        domains_added = len(list_domain_small_module_start)

        while count_module < numb_modules :
            if domain_start >= list_modules_few_domains_start[count_module] and domain_end <= list_modules_few_domains_end[count_module] :
                list_domain_small_module_start.append(domain_start)
                list_domain_small_module_end.append(domain_end)
                count_domain_end += 1
                count_module = numb_modules
                domains_added = len(list_domain_small_module_start)
            elif domain_start >= list_modules_few_domains_end[count_module] :
                count_module +=  1
            elif domain_start < list_modules_few_domains_start[count_module] :
                count_domain_end += 1
                count_module = numb_modules
            elif domain_start >= list_modules_few_domains_start[count_module] and domain_start <= list_modules_few_domains_end[count_module] and domain_end >= list_modules_few_domains_end[count_module]:
                domains_added = numb_domains_in_module
                count_module = numb_modules
            if domains_added == numb_domains_in_module :
                count_end_gap = 1
                for domain_end in list_domain_small_module_end :
                    count_numb_gap = 1
                    try :
                        gap_end = list_domain_small_module_start[count_end_gap] - 1
                    except :
                        continue
                    gap_start = domain_end
                    gap_id = str(gap_start)
                    # str_count_id_module = str(count_id_module)
                    # str_count_numb_gap = str(count_numb_gap)
                    # type_id = gap_id+".in_Module <= 3 domains"+str_count_id_module+"-"+str_count_numb_gap
                    type_id = gap_id+".in_Module <= 3 domains"

                    list_assign_domain_type.append(type_id)
                    list_all_gaps_start.append(gap_start)
                    list_all_gaps_end.append(gap_end)
                    count_numb_gap += 1
                    count_end_gap +=  1
                    
                    
    return list_all_gaps_start, list_all_gaps_end, list_assign_domain_type

def filter_gaps_in_one_cds(list_specific_gaps_start, list_specific_gaps_end, list_specific_gaps_type, list_cds_start, list_cds_end):
    count_end_gene = 0
    gaps_in_one_cds_start = list()
    gaps_in_one_cds_end = list()
    gaps_in_one_cds_type = list()
    for gene_start in list_specific_gaps_start :
        count_cds = 0
        numb_cds = int(len(list_cds_start))
        while count_cds < numb_cds :
            cds_start = list_cds_start[count_cds]
            cds_end = list_cds_end[count_cds]
            gene_end = list_specific_gaps_end[count_end_gene]
            gene_type = list_specific_gaps_type[count_end_gene]
            if cds_start <= gene_start <= cds_end and cds_start <= gene_end <= cds_end:
                gaps_in_one_cds_start.append(gene_start)
                gaps_in_one_cds_end.append(gene_end)
                gaps_in_one_cds_type.append(gene_type)
                count_cds = numb_cds
            else:
                count_cds +=1
        count_end_gene += 1
        
    return gaps_in_one_cds_start, gaps_in_one_cds_end, gaps_in_one_cds_type

def collect_gaps(list_specific_gaps_start, list_specific_gaps_end, list_specific_gap_type, list_collect_trans_gaps_start, list_collect_trans_gaps_end, list_collect_gap_types): 
    count_gap_end = 0
    for gap_start in list_specific_gaps_start:
        gap_end = list_specific_gaps_end[count_gap_end]
        gap_type = list_specific_gap_type[count_gap_end]
        list_collect_trans_gaps_start.append(gap_start)
        list_collect_trans_gaps_end.append(gap_end)
        list_collect_gap_types.append(gap_type)
        count_gap_end += 1
        
    return list_collect_trans_gaps_start, list_collect_trans_gaps_end, list_collect_gap_types

def sort_gap_types(list_sort_reference, list_to_get_sorted):
    list_to_get_sorted_done = list()
    numb_gaps = len(list_all_trans_gaps_start)
    for gap_start in list_sort_reference:
        count_gap = 0

        while count_gap < numb_gaps:
            gap_type_full = str(list_to_get_sorted[count_gap])
            gap_type_split = gap_type_full.split(".")
            gap_number = int(gap_type_split[0])
            if gap_start == gap_number:
                list_to_get_sorted_done.append(gap_type_full)
                list_to_get_sorted[count_gap] = 0
                count_gap = gap_number
            else:
                count_gap +=1
                
    return list_to_get_sorted_done

def get_long_gaps(list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gaps_type):
    list_long_trans_gaps_start = list()     #Liste mit Starts aller gaps die >600 bp sind
    list_long_trans_gaps_end = list()       #Liste mit Enden aller gaps die >600 bp sind
    list_long_gaps_type = list()
    count_gap = 0
    count_domains = len(list_all_trans_gaps_start)     #Anzahl der Lücken in Sequenz
    while count_gap < count_domains :
        pos_start_gap = list_all_trans_gaps_start[count_gap]
        pos_end_gap = list_all_trans_gaps_end[count_gap]
        type_gap = list_all_gaps_type[count_gap]
        length_gap = pos_end_gap - pos_start_gap
        if length_gap > 600:      #nur speichern seq wenn lücke größer als 0 ist
            list_long_trans_gaps_start.append(pos_start_gap)
            list_long_trans_gaps_end.append(pos_end_gap)
            list_long_gaps_type.append(type_gap)
        count_gap += 1
  
    return list_long_trans_gaps_start, list_long_trans_gaps_end, list_long_gaps_type

def delete_duplicate_gaps(list_gaps_start, list_gaps_end, list_gap_type):
    list_no_duplicate_gaps_start = list()
    list_no_duplicate_gaps_end = list()
    list_no_duplicate_gaps_type = list()
    numb_gaps = int(len(list_gaps_start))
    count_next_gap = 1
    count_gap = 0
    gap_test = 0
    while count_next_gap < numb_gaps :
        gap_start = list_gaps_start[count_gap]
        gap_end = list_gaps_end[count_gap]
        gap_type = list_gap_type[count_gap]
        gap_type_next = list_gap_type[count_next_gap]
        gap_start_next = list_gaps_start[count_next_gap]
        while gap_start == gap_start_next :
            gap_type = list_gap_type[count_gap]
            gap_type_next = list_gap_type[count_next_gap]
            gap_type = gap_type + " AND " + gap_type_next
            count_gap +=1
            count_next_gap +=1           
            if count_next_gap < numb_gaps:
                gap_start = list_gaps_start[count_gap]
                gap_end = list_gaps_end[count_gap]
                gap_start_next = list_gaps_start[count_next_gap]
            else:
                gap_start_next -=1
                gap_test = gap_start - 1            
        if gap_start_next == gap_test:
            count_gap -=1      
        gap_start = list_gaps_start[count_gap]
        gap_end = list_gaps_end[count_gap]
        # gap_type = list_gap_type[count_gap]
        list_no_duplicate_gaps_start+=[gap_start]
        list_no_duplicate_gaps_end+=[gap_end]
        list_no_duplicate_gaps_type+=[gap_type]
        count_gap +=1
        count_next_gap +=1
    if gap_test == 0:
        gap_start = list_gaps_start[count_gap]
        gap_end = list_gaps_end[count_gap]
        gap_type = list_gap_type[count_gap]
        list_no_duplicate_gaps_start+=[gap_start]
        list_no_duplicate_gaps_end+=[gap_end]
        list_no_duplicate_gaps_type+=[gap_type] 
        
    number_long_gaps = len(list_no_duplicate_gaps_start)
    
    return list_no_duplicate_gaps_start, list_no_duplicate_gaps_end, list_no_duplicate_gaps_type, number_long_gaps

def annotate_gene_seq_and_translation_gaps_save_seq_cds(list_long_trans_gaps_start, list_cds_in_region_start, list_cds_in_region_end, list_translation_cds): 
    count = 0
    count_domains = len(list_long_trans_gaps_start)
    #Speichern Sequenzen aus Lücken in Liste!!!
    list_gap_seq = list()
    while count < count_domains :
    #while count < 3 :
        pos_start_gap = int(list_long_trans_gaps_start[count] - 1)
        pos_end_gap = int(list_long_trans_gaps_end [count])
    #    print(repr(gbrecord.seq))
        #print(f"Seq Gap: {pos_start_gap}..{pos_end_gap}", repr(gbcontent.seq[pos_start_gap:pos_end_gap]))
        list_gap_seq.append(gbcontent.seq[pos_start_gap:pos_end_gap])
        count += 1
    
    list_aa_seq_gene = list()
    list_aa_seq_cds_with_gap = list()
    count_end_gene = 0
    for gene_start in list_long_trans_gaps_start :
        count_cds = 0
        numb_cds = int(len(list_cds_in_region_start))
        while count_cds < numb_cds :
            cds_start = list_cds_in_region_start[count_cds]
            cds_end = list_cds_in_region_end[count_cds]
            gene_end = list_long_trans_gaps_end[count_end_gene]
            if cds_start <= gene_start <= cds_end and cds_start <= gene_end <= cds_end:
                cds_domain = list_cds_features[count_cds]
                aa_seq_cds = list_translation_cds[count_cds]
                len_aa_seq_cds = len(aa_seq_cds)
                direction = cds_domain.strand 
                numb_aa_gene = int(((gene_end - gene_start) / 3) + 1)
                if direction > 0 :      #berechnen der Stelle an der AS Seq des neuen gens innerhalb der vorhanden CDS beginnt (forward strand)        
                    pos_first_aa = int(((gene_start - cds_start)/3)+1)
                    pos_last_aa = int(pos_first_aa + numb_aa_gene)
                    aa_seq_gene = aa_seq_cds[pos_first_aa:pos_last_aa]
                    aa_seq_len = len(aa_seq_gene)
                    list_aa_seq_gene.append(aa_seq_gene)
                    list_aa_seq_cds_with_gap.append(aa_seq_cds)
                    count_cds += 1
                else :
                    distance_cds_gene = int(((gene_start - cds_start)/3)-1)
                    pos_last_aa = int(len_aa_seq_cds - distance_cds_gene)
                    pos_first_aa = int(pos_last_aa - numb_aa_gene)
                    aa_seq_gene = aa_seq_cds[pos_first_aa:pos_last_aa]
                    aa_seq_len = len(aa_seq_gene)
                    list_aa_seq_gene.append(aa_seq_gene)
                    list_aa_seq_cds_with_gap.append(aa_seq_cds)
                    count_cds += 1

            
            else :
                count_cds += 1
                
        count_end_gene += 1
    
    return list_gap_seq, list_aa_seq_gene, list_aa_seq_cds_with_gap 

def generate_fasta_files(list_long_trans_gaps_start, list_long_trans_gaps_end, list_gap_seq, list_aa_seq_gene, list_aa_seq_cds_with_gap, id_file, list_gap_types_split):
    num_aa_seq = int(len(list_aa_seq_gene))
    num_gaps = len(list_long_trans_gaps_start)
    #id_file = gbcontent.id
    count_fasta = 1
    count = 0
    for create_fasta in list_gap_seq :
        if num_aa_seq > num_gaps :
            break
        description_gaps = "gap type:" , list_gap_types_split[count] ,"gap start:", list_long_trans_gaps_start[count], "gap end:", list_long_trans_gaps_end[count]
        str_description = str(description_gaps)
        gap_record = SeqRecord(create_fasta, id=id_file, description = str_description)
        str_count_fasta = str(count_fasta)
        gap_name = "\\Studium\\A_Master\\3. Semester\\trans-AT PKS\\GeneratedFiles\\Fasta\\"+id_file+"_Gene_Seq_Gap_"+str_count_fasta+".fasta"
        SeqIO.write(gap_record, gap_name, "fasta")
        
        description_cds = "gap type:" , list_gap_types_split[count] , "gap start:", list_long_trans_gaps_start [count], "gap end:", list_long_trans_gaps_end[count]
        str_cds_description = str(description_cds)
        aa_seq = list_aa_seq_gene[count]
        cds_record = SeqRecord(Seq(aa_seq), id=id_file, description=str_cds_description)
        cds_name = "\\Studium\\A_Master\\3. Semester\\trans-AT PKS\\GeneratedFiles\\Fasta\\"+id_file+"_AASeq_Gap_"+str_count_fasta+".fasta"
        SeqIO.write(cds_record, cds_name, "fasta")
        
        description_cds = "CDS Seq zur Lücke - gap type:" , list_gap_types_split[count] , "gap start:", list_long_trans_gaps_start [count], "gap end:", list_long_trans_gaps_end[count]
        str_cds_description = str(description_cds)
        aa_seq = list_aa_seq_cds_with_gap[count]
        whole_cds_record = SeqRecord(Seq(aa_seq), id=id_file, description=str_cds_description)
        whole_cds_name = "\\Studium\\A_Master\\3. Semester\\trans-AT PKS\\GeneratedFiles\\Fasta\\"+id_file+"_AASeq_CDS_of_the_Gap_"+str_count_fasta+".fasta"
        SeqIO.write(whole_cds_record, whole_cds_name, "fasta")
        count_fasta += 1
        count += 1

    return num_aa_seq, num_gaps

def generate_genbank_file(gbcontent, list_aa_seq_gene, list_cds_features, list_cds_in_region_start, list_cds_in_region_end, list_long_trans_gaps_start, list_long_trans_gaps_end):
    #Erstellen liste gaps_added wo zunächste alle features des ursprünglichen gbfiles drin sind
    gaps_added = list()
    for gbcontent in SeqIO.parse(gbfile, "genbank") :
                add_feature = gbcontent.features
    for feature in add_feature :   #loopen durch add_feature damit features die variablen haben, braucht 2. for loop :)
        gaps_added.append(feature)

    num_old_feats = len(gaps_added)

    #jetzt anhängen der Gaps als Gen feature und mit CDS feature des Gens an liste in der alle feature des gbk files schon drin sind
    count = 0
    count_cds = 0
   
    if num_aa_seq > num_gaps :
        count = num_gaps
    while count < num_gaps :
        try: 
            cds_feature = list_cds_features[count_cds]
        except:
            count = num_gaps
        if count < num_gaps :
            # c_start_cds = list_cds_in_region_start[count_cds]
            # c_end_cds = list_cds_in_region_end[count_cds]
            # c_start_gap = list_long_trans_gaps_start[count]
            if list_cds_in_region_start[count_cds] <= list_long_trans_gaps_start[count] <= list_cds_in_region_end[count_cds]:
                strand = cds_feature.strand
                if strand > 0 :
                    strand_feat = 1
                else:
                    strand_feat = -1
                try:
                    gene_feature = SeqFeature(FeatureLocation(list_long_trans_gaps_start[count], list_long_trans_gaps_end[count]), type ="gene", strand = strand_feat)
                except:
                    print("erstellen Seq feature gestoppt")
                    count = num_gaps
                gaps_added.append(gene_feature)
                aa_seq = list_aa_seq_gene[count]
                qualifier = {"translation":aa_seq}
                cds_feature = SeqFeature(FeatureLocation(list_long_trans_gaps_start[count], list_long_trans_gaps_end[count]), type ="cds", qualifiers=qualifier, strand = strand_feat)
                gaps_added.append(cds_feature)

                count += 1
            else :
                count_cds += 1
                
    num_added_feats = len(gaps_added)
    added_feats = int((num_added_feats - num_old_feats) /2)
   
    mol_type = {"molecule_type":"DNA"} #einführen Molekültyp, weil der für erstellen gbk file gebraucht wird
    #erstellen genbankfiles zur überprüfung der eigenen ergebnisse :)
    base_record = SeqRecord(gbcontent.seq, id = id_file, name = gbcontent.name, annotations=mol_type) #erstellen SeqRecord für gbk file
    for feats in gaps_added :
        base_record.features.append(feats)  #anhängen feature an base_record variable
    #Erstellen Datei: SeqIO.write(InhaltDokument, NameDokument, "TypDokument")
    output_name = "\\Studium\\A_Master\\3. Semester\\trans-AT PKS\\GeneratedFiles\\Genbank\\GenBank_"+id_file+"_addedgaps_"+".gbk"
    SeqIO.write(base_record, output_name, "genbank")   #erstellen gbkfile
    
    return added_feats

file_count = 0
count_files_with_gaps = 0
all_added_feats = 0       
for gbfile in files :
    file_count += 1
    anno, gbfeature, gbcontent, id_file = retrieve_features_from_genbankfile(gbfile)
    print(f"\nFile Nr: {file_count}: {gbfile}")
    
 #           Erstellen Liste mit allen CDS in Gensequenz  
    list_subregion_start, list_subregion_end, list_cds_in_region_start, list_cds_in_region_end, list_cds_features, list_all_domains_start, list_all_domains_end, list_domain_features, list_KS_zero_start, list_KS_zero_end, list_modules_few_domains_len, list_modules_few_domains_start, list_modules_few_domains_end, list_ACP_beta_start, list_ACP_beta_end, list_PP_start, list_PP_end = get_features(gbfeature)
 
    list_duplicate_domains_start, list_duplicate_domains_end, list_duplicate_type = generate_gaps_duplicate_domains(list_domain_features)

#       Speichern der Translation der CDS die in Subregions mit PKS/Polyketide sind
    list_translation_cds = get_translation_cds(gbfeature)
    
    #list_KS_zero_domain_end = cut_domains_longer_than_one_cds(list_cds_in_region_start, list_cds_in_region_end, list_KS_zero_start, list_KS_zero_end)
        
    list_KS_zero_gaps_start, list_KS_zero_gaps_end, list_KS_type = generate_gaps_before_and_after_domains(list_KS_zero_start, list_all_domains_start, list_all_domains_end, list_domain_features)
       
    list_ACP_beta_gaps_start, list_ACP_beta_gaps_end, list_ACP_type = generate_gaps_before_and_after_domains(list_ACP_beta_start, list_all_domains_start, list_all_domains_end, list_domain_features)

    list_PP_gaps_start, list_PP_gaps_end, list_PP_type = generate_gaps_before_and_after_domains(list_PP_start, list_all_domains_start, list_all_domains_end, list_domain_features)

    list_small_module_gaps_start, list_small_module_gaps_end, list_small_module_type = generate_small_module_gaps(list_modules_few_domains_len, list_modules_few_domains_start, list_modules_few_domains_end, list_all_domains_start, list_all_domains_end)

    list_small_module_gaps_start, list_small_module_gaps_end, list_small_module_type = filter_gaps_in_one_cds(list_small_module_gaps_start, list_small_module_gaps_end, list_small_module_type, list_cds_in_region_start, list_cds_in_region_end)
    
    list_collect_gaps_start = list()
    list_collect_gaps_end = list()
    list_collect_gap_types = list()        

    list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types = collect_gaps(list_KS_zero_gaps_start, list_KS_zero_gaps_end, list_KS_type, list_collect_gaps_start, list_collect_gaps_end, list_collect_gap_types)
    
    list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types = collect_gaps(list_ACP_beta_gaps_start, list_ACP_beta_gaps_end, list_ACP_type, list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types)
       
    list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types = collect_gaps(list_small_module_gaps_start, list_small_module_gaps_end, list_small_module_type, list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types)
   
    list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types = collect_gaps(list_duplicate_domains_start, list_duplicate_domains_end, list_duplicate_type, list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types)
  
    list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_trans_gaps_type = collect_gaps(list_PP_gaps_start, list_PP_gaps_end, list_PP_type, list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_gap_types)  
  
    list_all_trans_gaps_start.sort()
    list_all_trans_gaps_end.sort()
    
    list_all_trans_gaps_type = sort_gap_types(list_all_trans_gaps_start, list_all_trans_gaps_type)
            
    list_cut_gaps_start, list_cut_gaps_end = cut_domains_longer_than_one_cds(list_cds_in_region_start, list_cds_in_region_end, list_all_trans_gaps_start, list_all_trans_gaps_end, list_all_trans_gaps_type)
    
    list_long_trans_gaps_start, list_long_trans_gaps_end, list_long_gaps_type = get_long_gaps(list_cut_gaps_start, list_cut_gaps_end, list_all_trans_gaps_type)
    
    try_list = list()
    try :
        try_list.append(list_long_trans_gaps_start[0])
    except:
        print("Keine Gaps >600pb")
        continue
    count_files_with_gaps +=1
    
    list_long_trans_gaps_start, list_long_trans_gaps_end, list_long_gaps_type, number_long_gaps = delete_duplicate_gaps(list_long_trans_gaps_start, list_long_trans_gaps_end, list_long_gaps_type)    
    
    list_gap_types_split = list()
    list_multiple_types = list()
    #list_long_gaps_type = ["123.abc","1222.abc AND 456.lol","890.hihihi AND 1222.abc AND 456.lol"]
    #Sortieren der Typen der Lücken!
    for gap_type_full in list_long_gaps_type:
        if "AND" in gap_type_full:
            list_multiple_types = list()
            gap_type_split = gap_type_full.split("AND")
            numb_types = len(gap_type_split)
            multiple_types = " "
            for types in gap_type_split:
                split_type = types.split(".")
                type_name = split_type[1]
                list_multiple_types.append(type_name)
            count_type = 0
            count_type_next = 1
            combined_type = " "
            while count_type_next < numb_types:
                combined_type_now = list_multiple_types[count_type]
                combined_type_next = list_multiple_types[count_type_next]
                if "AND" in combined_type:                    
                    combined_type = combined_type + "AND " + combined_type_next
                else:
                    combined_type = combined_type_now + "AND " + combined_type_next
                count_type +=1
                count_type_next +=1    
            list_gap_types_split.append(combined_type)
        else:
            gap_type_split = gap_type_full.split(".")
            gap_type = gap_type_split[1]
            list_gap_types_split.append(gap_type)
    
   
    
    list_gap_seq, list_aa_seq_gene, list_aa_seq_cds_with_gap = annotate_gene_seq_and_translation_gaps_save_seq_cds(list_long_trans_gaps_start, list_cds_in_region_start, list_cds_in_region_end, list_translation_cds)

    #Exportieren jeder Lücke als Fasta Files mit ID und Beschreibung, ein FIle mit gen seq und ein file mit AS seq
    #Erstellen Variable die in den Namen der neu erstellten datei kommen soll --> wiedererkennung datei

    num_aa_seq, num_gaps = generate_fasta_files(list_long_trans_gaps_start, list_long_trans_gaps_end, list_gap_seq, list_aa_seq_gene, list_aa_seq_cds_with_gap, id_file, list_gap_types_split)

    added_feats = generate_genbank_file(gbcontent, list_aa_seq_gene, list_cds_features, list_cds_in_region_start, list_cds_in_region_end, list_long_trans_gaps_start, list_long_trans_gaps_end)
    
    all_added_feats = all_added_feats + added_feats
    
    print(f"In {id_file} sind {number_long_gaps} Gaps >600 bp, es wurden {added_feats} Gaps an das Genbank File angehängt")
    
    #Feature SeqIO.write(InhaltDokument, NameDokument, "TypDokument")
    #Importieren der Gaps als "gene" in genbankfile
print("\n        *********")
print(f"\nEs wurden insgesamt {all_added_feats} Gaps mit einer Länge von >600bp in {count_files_with_gaps} von {file_count} Files generiert.")

#codeword = input("Codeword? ")
#if "music" in codeword or "musik" in codeword :
#    print(music)
#elif "bio" in codeword or "Bio" in codeword:
#    print(bio)
#elif "trans" in codeword or "Trans" in codeword:
#    print(trans_PKS)
#else :
#    print("Sorry, u r wrong.")
