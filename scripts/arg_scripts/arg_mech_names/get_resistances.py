import sys
import os

folder='CDC2'
argpid='99'


obofile ='/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/aro.obo.agk'

def get_stanzas(obofile):
	stanzas = []
	with open(obofile) as obo:
		stanza = ''
		for line in obo:
			if line != '\n':
				stanza = stanza + line
			if line == '\n':
				stanzas.append(stanza)
				stanza = ''
	return stanzas


arodict = {}
aroisadict = {}
stanzas = get_stanzas(obofile)
outhandle = 'aro_name_drug.txt'
with open(outhandle,'w') as outfile:
	for stanza in stanzas:
		aroids = []
		if '[Term]' in stanza:
			aroid = stanza.split('\nid: ')[1].split('\n')[0]
			aroids.append(aroid)
			aroname = stanza.split('\nname: ')[1].split('\n')[0]
			if 'is_a' in stanza:
				aroisas = stanza.split('\nis_a: ')
				for links in aroisas[1:]:
					aroisa = links.split(' ! ')[0]
					aroids.append(aroisa)
			if 'relationship' in stanza:
				arorels = stanza.split('\nrelationship: ')
				for links in arorels[1:]:
					relationship = links.split(' ARO')[0]
					if relationship == 'confers_resistance_to_drug':
						drug = links.split(' ! ')[1].split('\n')[0]
						for aroid in aroids:
							outfile.write('{0}\t{1}\t{2}\n'.format(aroid,aroname,drug))


