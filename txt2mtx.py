import sys
import os
import porter

data_file = sys.argv[1]

do_stemming = False
stemmer = porter.PorterStemmer()
stemmer_mapping = {}

def truncate_word(word):
	start = 0
	while start < len(word) and word[start].isalnum() == False:
		start += 1
	end = len(word)
	while end > start and word[end-1].isalnum() == False:
		end -= 1
	truncated = word[start:end].lower()
	for letter in truncated:
		if letter.isalpha():
			break
	else:
		return ''
	try:
		truncated.decode('ascii')
	except UnicodeDecodeError:
		return ''
	if do_stemming == True:
		if len(truncated) == 0:
			return ''
		else:
			stemmed = stemmer.stem(truncated, 0, len(truncated)-1)
			stemmer_mapping[truncated] = stemmed
			return stemmed
	else:
		return truncated

def read_txt(line):
	#line = line[0:-1].split('\t')[-1]
	bag_words = dict()
	row = line.split()
	for word in row:
		truncated = truncate_word(word)
		if truncated != '':
			if truncated in bag_words:
				bag_words[truncated] += 1
			else:
				bag_words[truncated] = 1
	return bag_words

stop_list = set()
f_stop = open('english.stop')
for line in f_stop:
	word = line[0:-1]
	if do_stemming == True:
		word = stemmer.stem(word, 0, len(word)-1)
	stop_list.add(word)
f_stop.close()

bag_words = dict()
f_input = open(data_file + '.txt')
#doc_file = open('documents.txt', 'w')
doc_count = 0
for line in f_input:
	bag_words_one = read_txt(line)
#	doc_file.write(filename + '\n')
	for word in bag_words_one:
		if word in bag_words:
			bag_words[word] += bag_words_one[word]
		else:
			bag_words[word] = bag_words_one[word]
	doc_count += 1
	if doc_count % 10000 == 0:
		print doc_count, 'documents processed...'
#doc_file.close()
f_input.close()

voc_file = open('vocabulary.txt', 'w')
word_map = dict()
count = 0
for word in sorted(bag_words):
	if word in stop_list or bag_words[word] < 5:
		continue
	voc_file.write(word + '\t' + str(bag_words[word]) + '\n')
	word_map[word] = count
	count += 1
voc_file.close()

f_mtx = open(data_file + '.mtx', 'w')
f_mtx.write('%%MatrixMarket matrix coordinate real general\n%\n')
f_mtx.write(str(len(word_map)) + ' ' + str(doc_count) + ' 0\n')
doc_count = 0
f_input = open(data_file + '.txt')
for line in f_input:
        bag_words_one = read_txt(line)
        for word in sorted(bag_words_one):
                try:
                        word_idx = word_map[word]
                        f_mtx.write(str(word_idx+1) + ' ' + str(doc_count+1) + ' ' + str(bag_words_one[word]) + '\n')
                except KeyError:
                        continue
        doc_count += 1
        if doc_count % 10000 == 0:
                print doc_count, 'documents processed ......'
f_input.close()
f_mtx.close()

f_stem = open('stemmer_mapping.txt', 'w')
for word in sorted(stemmer_mapping.iteritems(), key=lambda i: i[1]):
	f_stem.write(word[1] + '\t' + word[0] + '\n')
f_stem.close()
