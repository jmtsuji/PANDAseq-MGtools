#Enter INPUTFILE, OVERLAP, and BUFFER values at lines 54-56
"""Copyright Katherine Naismith, 2016
**Still an EARLY TESTING VERSION -- use at your own risk.
Contact K. Naismith if you have usage questions -- kenaismith@uwaterloo.ca
    
This is designed to read in a hit table text
file output by the online database.

It is assumed that results with the same label 
occur on consecutive lines within the same query 

-------------OUTPUT----------------
Match Found: some specimen (label given) matched the input sequence
	with bases matching on BOTH sides of the stitched part of the sequence
	
Overlap: some 'left side' and 'right side' match on the same specimen
	overlap, but no match spans the overlapping portion. The overlap 
	value closest to the OVERLAP parameter given is returned, along with
	a label on which it is attained.

Gap: the closest correctly oriented 'left side' and 
	'right side' matches of the sequence across all specimens had a 
	gap between the end of the 'left side' and the beginning 
	of the 'right side'. The smallest gap and the specimen on which it
	was attained are given. A gap value of 0 means the 'left side'
	and 'right side' land end-to-end with no overlap.
	
	i.e.	
	
		left start ---left end
							  |-gap-|
									 right start------right end	
	
Results Not Comparable: There was no pair of 'left side' and 'right side'
	matches, on any specimen returned by the query, where the 'left side'
	and 'right side' had the same orientation and where the 'left side' 
	started before the 'right side'
	i.e. All pairs of left and right matches for the query had
	
		right start ----------right end
						left end ------------left start
	
		OR
	
		right start --------------right end 
					left start here--------------
					or further --->		
"""

#you don't need the next line, since all functions are in here
#from chimeralibrary import *

#Set the input file name, and overlap size here
INPUTFILE = 'RID-Alignment.txt'
OVERLAP = 10
BUFFER = 10



#FUNCTIONS--------------------------------------------------------
'''getchar takes a file as input
and returns the next character in the file.
'1' is the number of bytes to be read
by f.read(size).
'''
def getchar(INPUTFILE):
	char = INPUTFILE.read(1)
	return char 
	
'''getblock takes a file as input
and returns the block of text starting with the
current character and ending with the character 
before the next tab character. Note that
getblock still reads the next tab character.'''
def getblock(INPUTFILE):
	block = ''
	char = getchar(INPUTFILE)
	while (char != '\t'):
		block = block + char
		char = getchar(INPUTFILE)
	return(block)

'''handleline takes as input: the input and output files, the current
left and right match lists and the start & end buffer values. 
It assumes the last character read on the line is the tab before the sample start value.
handleline updates our lists of right and left matches
to include the info on this line
numbers in the backward list are negative
min_tag = True if no match is found that spans the middle of the stitched segment.
min_tag = False if sugh a match is found (i.e. no min value to be generated)'''
def handleline(INPUTFILE, OUTPUTFILE, BUFFER, label, L_forward_list, L_backward_list, R_forward_list, R_backward_list, L_START, L_END, R_START, R_END, OVERLAP):
	#print 'handle line'
	min_tag = True
	#Get query start and end. q_start is ALWAYS less than q_end, by choice of file type.
	q_start = int(getblock(INPUTFILE))
	q_end = int(getblock(INPUTFILE))
		#get sample start and end
	s_start = int(getblock(INPUTFILE))
	s_end = int(getblock(INPUTFILE))
#Add start and end values to appropriate list
	if (L_START <= q_start <= L_END - BUFFER):
		if (R_START + BUFFER <= q_end <= R_END): #match starts in left and ends in right
			#OUTPUTFILE.write('Alignment Succeeded\t%s \n' %label)
			min_tag = False
			#OUTPUTFILE.write('q_start = %d, q_end = %d \n' %(q_start, q_end))
			INPUTFILE.readline()
			flag = getchar(INPUTFILE)
			while (flag!= '#' and flag > 32 ):
				INPUTFILE.readline()
				flag = getchar(INPUTFILE)
				#print 'ord(flag) = %d' %ord(flag)
				#print'loop'
		else:#match starts in left and ends in left
			if (s_start < s_end):
				L_forward_list.append(s_start - (q_start - L_START))
				L_forward_list.append(s_end + (L_END + OVERLAP - q_end))
			else:
				L_backward_list.append((-1)*(s_start + q_start - L_START))
				L_backward_list.append((-1)*(s_end - (L_END + OVERLAP - q_end)))
			INPUTFILE.readline()
			flag = getchar(INPUTFILE)

	elif (R_START + BUFFER <= q_end <= R_END): #then the match starts in the right piece, and stays there
		if (s_start < s_end):
			R_forward_list.append(s_start - (q_start + OVERLAP - R_START))
			R_forward_list.append(s_end + (R_END - q_end))
		else:
			R_backward_list.append((-1)*(s_start + q_start + OVERLAP - R_START))
			R_backward_list.append((-1)*(s_end - (R_END - q_end)))
		INPUTFILE.readline()
		flag = getchar(INPUTFILE)
	else: #match starts in buffer/overlap
		if (q_end <= L_END - BUFFER): #in left backward list
			L_backward_list.append((-1)*(s_start + q_start - L_START))
			L_backward_list.append((-1)*(s_end - (L_END + OVERLAP - q_end)))
		elif (q_end >= R_START + BUFFER): #in right forward list
			R_forward_list.append(s_start - (q_start + OVERLAP - R_START))
			R_forward_list.append(s_end + (R_END - q_end))		
		INPUTFILE.readline()
		flag = getchar(INPUTFILE)
	return (flag, min_tag)

'''mingap takes a left list and a right list, each with even values "start" 
values and odd values "end" values. These lists must have "start" values 
less than "end" values for each pair in positions 2i, 2i+1. 
mingap returns the diff/overlap distance between
a left "end" value and a right "start" value that is closest to 1.
The distance is given as negative if there is overlap, and positive 
if there is a diff. In the case of a tie between values, an arbitrary
value is given'''
def mingap(L_list, R_list, OVERLAP):	
	i=0
	diff = 40000000000
	while (2*i < len(L_list)):
		j=0
		while(2*j < len(R_list)):
			if (L_list[2*i] < R_list[2*j]):
				tmp_diff = R_list[2*j] - L_list[2*i+1] - 1 
				diff = compare(tmp_diff, diff, OVERLAP)	
			j=j+1
		i = i+1
	return(diff)
	
'''compare(tpm_diff, diff, OVERLAP) compares the possible new difference value (tmp_diff) to 
the previous difference value (diff). 
The higher a value ranks in the following list, the better it is:
1) negative value as close to -OVERLAP as possible 
2) positive value closest to zero'''
def compare(tmp_diff, diff, OVERLAP):
	if (diff < 0): #tmp_diff represents an overlap
		if (tmp_diff < 0): #diff also represents an overlap
			if (abs(tmp_diff - (OVERLAP - 1)) < abs(tmp_diff - (OVERLAP - 1))):#if tmp_diff is a 'better' overlap
				diff = tmp_diff
		#else: tmp_diff is a gap. do nothing
	else: #diff represents a gap 
		if (tmp_diff <0): #tmp_diff is an overlap
			diff = tmp_diff 
		else: #tmp_diff is also a gap
			if(tmp_diff < diff): #if tmp_diff is a smaller gap
				diff = tmp_diff
	return diff
			
''' start running handle_hashtag_lines when the most recent character
read is the first hashtag of the first line of a block of hashtags.
The function gets the query name, reads through all hashtagged lines
for the query, and ends treatment of the query if there are no results
'''			
def handle_hashtag_lines(INPUTFILE, OUTPUTFILE, char):
	hashtag = 1
	overcount = 0
	while (char == '#'):
		if (hashtag == 3):
			if overcount ==1:
				OUTPUTFILE.write('No Results\tNA\tNA\n')
				overcount = 0
			getchar(INPUTFILE)
			query_name = ''
			nxt = getchar(INPUTFILE)
			while nxt != '\n':
				query_name = query_name + nxt
				nxt = getchar(INPUTFILE)
			OUTPUTFILE.write('%s \t' %query_name)
			hashtag = 4
			char = getchar(INPUTFILE)
		elif (hashtag == 5):
			overcount = 1
			INPUTFILE.readline()
			char = getchar(INPUTFILE)
			if (ord(char) > 32):
				hashtag = 1
			else:#end of file
				OUTPUTFILE.write('No Results\tNA\tNA\n')
		else:
			INPUTFILE.readline()
			char = getchar(INPUTFILE)
			hashtag = hashtag + 1
	return(char)
	
'''	finishlabel gets the best diff (gap or overlap) value for the specimen,
and the new value if it is better than the previous diff for the query.'''
def finishlabel(L_forward_list, L_backward_list, R_forward_list, R_backward_list, query_diff, OVERLAP):
	forward_diff = mingap(L_forward_list, R_forward_list, OVERLAP)
	backward_diff = mingap(L_backward_list, R_backward_list, OVERLAP)
	label_diff = compare(forward_diff, backward_diff, OVERLAP)
	query_diff = compare(query_diff, label_diff, OVERLAP)
	return (query_diff)



#MAIN----------------------------------------------------------------

#names the output file
OUTPUTFILE = 'output_'+ INPUTFILE

#Set boundaries for left and right reads
L_START = 0
L_END = 200 - OVERLAP
R_START = 201
R_END = 401 - OVERLAP

#Number of tabs before 'important' types of data on a line of the hit table
LABEL_LOCATION = 1
SAMPLE_START_LOCATION = 6

#open file
chimera_input = open( INPUTFILE, 'r' )
chimera_output = open( OUTPUTFILE, 'w')


#Print header in output file
chimera_output.write('#Query\tResult\tSize(bases)\tSubject\n')

#First iteration of loop
char = getchar(chimera_input)

while (ord(char)>32): #char not white space
	#Do hashtag lines until we run out
	char = handle_hashtag_lines(chimera_input, chimera_output, char)
	if (ord(char) > 32): #char not white space
	# Now we have a (real) line that doesn't start with a hashtag. Get sample name (label).
		tabcount = 0
		while (tabcount < LABEL_LOCATION):
			getblock(chimera_input)
			tabcount = tabcount + 1
		label = getblock(chimera_input)
		tabcount = tabcount+1
	#go to sample_start
		while(tabcount < SAMPLE_START_LOCATION):
			getblock(chimera_input)
			tabcount = tabcount + 1
	#initialize start/end values, difference, lists.
		left_match_start = 0
		left_match_end = 0
		right_match_start = 40000000000
		right_match_end = 40000000000
	#initialize query_diff
		query_diff = right_match_start - left_match_end
	#the lists have even indices "start" values, odd indices "end" values
	#values in backward lists must be multiplied by -1 before storing.
		L_forward_list = []
		L_backward_list = []
		R_forward_list = []
		R_backward_list = []
	#deal with the line, add to list
		char, min_tag = handleline(chimera_input, chimera_output, BUFFER, label, L_forward_list, L_backward_list, R_forward_list, R_backward_list, L_START, L_END, R_START, R_END, OVERLAP)
		while ((ord(char) > 32 or ord(char)== 9) and char != '#'): 
			tabcount = 0
			getblock(chimera_input)
			tabcount = tabcount + 1
			newlabel = getblock(chimera_input)
			tabcount = tabcount + 1
			if (newlabel!=label): #new sample. get the label_min from the label we just finished. Then re-initialize label and lists
				query_diff = finishlabel(L_forward_list, L_backward_list, R_forward_list, R_backward_list, query_diff, OVERLAP)
				query_label = label
			#reinitialize label and lists
				label = newlabel
				while(tabcount < SAMPLE_START_LOCATION):
					getblock(chimera_input)
					tabcount = tabcount + 1

				left_match_start = 0
				left_match_end = 0
				right_match_start = 40000000000
				right_match_end = 40000000000
			

				#the lists have even indices "start" values, odd indices "end" values
				#values in backward lists must be multiplied by -1 before storing.
				L_forward_list = []
				L_backward_list = []
				R_forward_list = []
				R_backward_list = []
	
				char, min_tag = handleline(chimera_input, chimera_output, BUFFER, label, L_forward_list, L_backward_list, R_forward_list, R_backward_list, L_START, L_END, R_START, R_END, OVERLAP)
				
			else: #then we keep appending to our current lists
				
				while(tabcount < SAMPLE_START_LOCATION):
					getblock(chimera_input)
					tabcount = tabcount + 1
					
				char, min_tag = handleline(chimera_input, chimera_output, BUFFER, label, L_forward_list, L_backward_list, R_forward_list, R_backward_list, L_START, L_END, R_START, R_END, OVERLAP)

		if (char == '#' or ord(char) <= 32):
			query_diff = finishlabel(L_forward_list, L_backward_list, R_forward_list, R_backward_list, query_diff, OVERLAP)
			query_label = label
			if (min_tag ==False):
				chimera_output.write('Alignment Succeeded\tNA\t%s\n' %label)
			elif (query_diff < 40000000000):
				if (query_diff <0):
					chimera_output.write("Overlap\t%d\t%s\n" % (abs(query_diff), query_label))
				else:
					chimera_output.write("Gap\t%d\t%s\n" %(query_diff, query_label))
			else:
				chimera_output.write("Results Not Comparable\tNA\tNA\n")
			
print 'done with file %s' %INPUTFILE
