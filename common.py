#import Levenshtein 

class arithmetic():  
      
    def __init__(self):  
        pass  

    def levenshtein(self,first,second):  
        if len(first) > len(second):  
            first,second = second,first  
        if len(first) == 0:  
            return len(second)  
        if len(second) == 0:  
            return len(first)  
        first_length = len(first) + 1  
        second_length = len(second) + 1  
        distance_matrix = [range(second_length) for x in range(first_length)]   
        #print distance_matrix  
        for i in range(1,first_length):  
            for j in range(1,second_length):  
                deletion = distance_matrix[i-1][j] + 1  
                insertion = distance_matrix[i][j-1] + 1  
                substitution = distance_matrix[i-1][j-1]  
                if first[i-1] != second[j-1]:  
                    substitution += 1  
                distance_matrix[i][j] = min(insertion,deletion,substitution)  
        #print distance_matrix  
        return distance_matrix[first_length-1][second_length-1]
    
    def normal_leven(self, str1, str2):
        len_str1 = len(str1) + 1
        len_str2 = len(str2) + 1
        #create matrix
        matrix = [0 for n in range(len_str1 * len_str2)]
        #init x axis
        for i in range(len_str1):
            matrix[i] = i
        #init y axis
        for j in range(0, len(matrix), len_str1):
            if j % len_str1 == 0:
                matrix[j] = j // len_str1
        
        for i in range(1, len_str1):
            for j in range(1, len_str2):
                if str1[i-1] == str2[j-1]:
                    cost = 0
                else:
                    cost = 1
                matrix[j*len_str1+i] = min(matrix[(j-1)*len_str1+i]+1,
                          matrix[j*len_str1+(i-1)]+1,
                          matrix[(j-1)*len_str1+(i-1)] + cost)
        return matrix[-1]
    
    def checkErr(self,seq_s,seq_q,errorN=1):
           
        ##temp_N = Levenshtein.distance(seq_s,seq_q) 
        temp_N = self.levenshtein(seq_s,seq_q) ## Use the local function, in case the Levenshtein Package is not installed.
        if temp_N > errorN:
            return 0  ##failed, Error Num > limited N
        else:
            return 1 


    
    def checkError_hrf(self,seq_s,seq_q,errorN=1):
        seq_s_list = list(seq_s)
        seq_q_list = list(seq_q)
        s_l = len(seq_s_list)-1
        q_l = len(seq_q_list)-1
        error_n = 0
        q_pos = 0
        s_pos = 0
        while s_pos < s_l and q_pos < q_l:
            if error_n > errorN:
                break
            if seq_s_list[s_pos] != seq_q_list[q_pos]:
                error_n += 1
                if seq_s_list[s_pos+1] == seq_q_list[q_pos+1]:
                    ## mismatch
                    s_pos += 1
                    q_pos += 1
                    continue
    
                if seq_s_list[s_pos] == seq_q_list[q_pos+1]:
                    s_pos += 1
                    q_pos += 2
                    continue
                
                if seq_s_list[s_pos+1] == seq_q_list[q_pos]:
                    s_pos += 2
                    q_pos += 1
                    continue
            
            s_pos += 1
            q_pos += 1
            
        return 0 if error_n > errorN else q_pos
 
    