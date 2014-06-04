'''
Base classes+functions for hierarchy classification project

methods:
  training:
    build tree for determining BFS path
    for each node with children, compile average terms down both paths
  testing:
    follow decision tree, choosing best path at each node
    compile a list of nodes visited
    once end of tree is reached, check against all nodes visited to decide
    which node(s) the test doc belongs to.
'''

import simplejson as json
from math import sqrt

def first_load_hierarchy():
  # Run this after making any major changes. It compiles helper files to be
  #   used for easier analysis.
  print "Compiling hierarchy..."
  compile_hierarchy()
  print "Reading hierarchy from file..."
  nodes = read_hierarchy()
  print "Determining parent nodes..."
  compile_hierarchy_parents(nodes)
  print "Load Complete"
  return

def compile_hierarchy(rfname='hierarchy.txt', wfname='hierarchy_json.txt'):
  # Reads the hierarchy file and processes it into a JSON file for easy
  #   future reference
  nodes = {}
  with open(rfname) as f:
    while True:
      l = f.readline()
      if l=='':
        break
      l = l.split()
      if l[0] in nodes:
        nodes[l[0]].append(l[1])
      else:
        nodes[l[0]] = [l[1]]
  with open(wfname,'w') as wf:
    json.dump(nodes, wf)
  return

def read_hierarchy(rfname='hierarchy_json.txt'):
  # Reads the JSON hierarchy file, returns a directory of nodes
  nodes = {}
  with open(rfname) as f:
    nodes = json.loads(f.read())
  return nodes

def hierarchy_info(nodes):
  # Prints some info about the hierarchy. Edit as needed for testing.
  histocount = {}
  for n in nodes:
    c = len(nodes[n])
    if c in histocount:
      histocount[c] += 1
    else:
      histocount[c] = 1
  for h in histocount:
    print '%s nodes with exactly %s children' % (histocount[h], h)
  return

def compile_hierarchy_parents(nodes):
  # Determines which nodes have no parents. Writes list to 'parent_nodes.txt'.
  seen = set([])
  parents = set([])
  for n in nodes:
    if n not in seen:
      parents.add(n)
    for c in nodes[n]:
      seen.add(c)
      parents.discard(c)
  parents = list(parents)
  with open('parent_nodes.txt','w') as pf:
    json.dump(parents, pf)
  return

def load_parents(fname='parent_nodes.txt'):
  # Loads from the parent nodes file into a list variable.
  parents = []
  with open(fname) as f:
    parents = json.loads(f.read())
  return parents

def all_children(nodes, parent):
  # Returns a set of all children and sub-children of node 'parent'
  #   using a BFS algorithm
  children = set([])
  queue = [parent]
  while True:
    n = queue[0]
    for c in nodes[n]:
      if c not in children:
        queue.append(c)
    if len(queue) == 1:
      break
    queue = queue[1:]
  return list(children)

def compile_training_records(terms_list=None, fname='train-remapped.csv'):
  # Reads the training file
  # Creates a dictionary of simple categories:
  #   key: training record id, value: list of its categories.
  # Also creates a dict of terms and documents:
  #   key: term, value: list of document ids
  # Note: terms dict currently considers only existence of terms, not frequency
  # term_list allows for limiting the number of terms loaded into the dict.
  #   use a list from load_terms_counts with a numbered cutoff if desired.
  cats_dict = {}
  terms_dict = {}
  i = 10000001
  with open(fname) as rf:
    rf.readline() # clear header
    while True:
      l = rf.readline()
      if l=='':
        break
      nodes, terms = train_line_parser(l)
      cats_dict[str(i)] = nodes
      if terms_list == None:
        for t in terms:
          if t in terms_dict:
            terms_dict[t].append(str(i))
          else:
            terms_dict[t] = [str(i)]
      else:
        for t in terms:
          if t in terms_list:
            if t in terms_dict:
              terms_dict[t].append(str(i))
            else:
              terms_dict[t] = [str(i)]
      i += 1
  with open('training_categories.txt','w') as wfc:
    json.dump(cats_dict, wfc)
  with open('terms_to_records.txt','w') as wft:
    json.dump(terms_dict, wft)
  return

def load_training_records(fname='training_categories.txt'):
  # Loads the JSON dict of training records and their categories
  records = {}
  with open(fname) as f:
    records = json.loads(f.read())
  return records

def compile_category_reference_frequencies(records):
  # Creates a dict. key: category, value: number of references in training data
  #   Used to weight scores in the analyze() algorithm
  ref_freq = {}
  for r in records:
    for c in records[r]:
      if c in ref_freq:
        ref_freq[c] += 1
      else:
        ref_freq[c] = 1
  return ref_freq

def load_terms_records(fname='terms_to_records.txt'):
  terms = {}
  with open(fname) as f:
    terms = json.loads(f.read())
  return terms

def train_line_parser(line):
  # Takes a line of training data, processes it, and returns the useful info
  #   it contains: a list of nodes it belongs to and a dict of terms
  nodes = []
  terms = {}
  l = line.split()
  i = 0
  while ':' not in l[i]:
    # Extract nodes
    nodes.append(l[i].strip(','))
    i += 1
  for j in range(i,len(l)):
    t = l[j].split(':')
    terms[t[0]] = int(t[1])
  return nodes, terms

def test_line_parser(line):
  # Same as train_line_parser but for test data.
  terms = {}
  l = line.split()
  recordid = l[0].split(',')
  recordid = recordid[0]
  for i in range(1,len(l)):
    t = l[i].split(':')
    terms[t[0]] = int(t[1])
  return recordid, terms

def sum_training_terms(n=1000, tfile='train-remapped.csv'):
  # Reads the terms from the first n entries in the training data
  # Creates a master dict of all terms and their total counts
  # Will be useful for determining the relative rarity of a term.
  d = {}
  with open(tfile) as rf:
    rf.readline() # clear header
    if n=='all':
      while True:
        l = rf.readline()
        if l=='':
          break
        nodes, terms = train_line_parser(l)
        for t in terms:
          if t in d:
            d[t] += terms[t]
          else:
            d[t] = terms[t]
    else:
      for i in range(n):
        nodes, terms = train_line_parser(rf.readline())
        for t in terms:
          if t in d:
            d[t] += terms[t]
          else:
            d[t] = terms[t]
  with open('terms_counts.txt','w') as wf:
    json.dump(d,wf)
  return

def load_terms_counts(cutoff=None, fname='terms_counts.txt'):
  # Loads the terms_counts file into a dict.
  #   key: term, value: count (in all training records)
  termcount = {}
  with open(fname) as f:
    termcount = json.loads(f.read())
  if cutoff == None:
    return termcount
  else:
    limited_terms = {}
    for t in termcount:
      if termcount[t] <= cutoff:
        limited_terms[t] = termcount[t]
    return limited_terms

def terms_by_count(terms):
  # Takes a dict of terms (most likely from sum_training_terms) and determines
  #   the frequency of each count, creating a dict as follows:
  #   key: count, value: number of terms with that count
  counts = {}
  for t in terms:
    if str(terms[t]) in counts:
      counts[str(terms[t])] += 1
    else:
      counts[str(terms[t])] = 1
  return counts

def recompile_terms_cutoff(cutoff):
  # Reruns functions to compile terms lists+dicts with a cutoff threshold of
  #   term count. Run this before running analyze() to maximize efficiency.
  print 'Loading term counts...'
  terms = load_terms_counts(cutoff)
  print 'Compiling training records...'
  compile_training_records(terms)
  print 'Recompile complete!'
  return

def clean_prediction_file(fname='djsensei-predictions.csv'):
  with open(fname) as rf:
    with open('djsensei-predictions-2.csv','w') as wf:
      header = rf.readline()
      wf.write(header)
      while True:
        l = rf.readline()
        if l=='':
          break
        l = l.strip(' \n') + '\n'
        wf.write(l)
  return

def analyze(term_cutoff = 5000, \
            testing_count = 100, \
            fname='test-remapped.csv', \
            outname='djsensei-predictions.txt'):
  # The test analyzer! Loads dicts to work with, then reads in test data.

  # Load dicts for easy reference
  print 'Loading training records and categories...'
  rec_cat = load_training_records()
  print 'Loading category reference frequencies...'
  ref_freq = compile_category_reference_frequencies(rec_cat)
  print 'Loading terms and records...'
  term_rec = load_terms_records()
  print 'Loading terms and counts...'
  term_count = load_terms_counts(term_cutoff)

  # Read in test line
  print 'Reading test records...'
  with open(fname) as rf:
    rf.readline() # clear header
    with open(outname, 'w') as wf:
      wf.write('Id,Predicted\n')
      borked = [] # problem lines in the test record, come back to those later.
      i = 1
      while True: # for i in range(1,testing_count):
        l = rf.readline()
        if l == '':
          break
        if i%1000 == 0:
          print 'Reading test record %s...' % i

        pos_cat = {} # key: a possible category, value: a weighted factor
        recordid, terms = test_line_parser(l)

        # Analyze test line against existing records
        # Compile possible categories from matching terms
        for t in terms:
          if t in term_count:
            tc = term_count[t]
            for r in term_rec[t]:
              for c in rec_cat[r]:
                score = float(terms[t])/(sqrt(tc)*sqrt(ref_freq[c]))
                if c in pos_cat:
                  pos_cat[c] += score
                else:
                  pos_cat[c] = score

        # Sort possible categories by highest scores.
        tops = sorted(pos_cat, key = lambda k: pos_cat[k])
        tops.reverse()

        # Threshold for inclusion: greater than .5 * the best score
        last = 1
        while tops != []:
          if last == len(tops):
            break
          if pos_cat[tops[last]] / pos_cat[tops[0]] < .75 or last > 4:
            break
          else:
            last += 1
        tops = tops[:last]
        if tops == []:
          borked.append(recordid)

        # Write predicted categories to output file
        write_string = recordid + ','
        for t in tops:
          write_string += t + ' '
        write_string.strip(' ')
        wf.write(write_string+'\n')
        i += 1
  with open('borked_tests.txt','w') as bff:
    json.dump(borked,bff)
  return

if __name__ == "__main__":
  analyze()
