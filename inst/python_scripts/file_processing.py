# Copyright 2012 Oliver Serang
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

import sys
from lxml import etree

def get_assayID_sampleID_mass_height_area_from_XML(filename):
    xml_tree = etree.parse(filename)

    records = xml_tree.findall('.//record')
    for r in records:
        peakType = r.attrib.get('peakType')
        
        if peakType != 'A':
            continue

        assayID = r.attrib.get('assayId')
        sampleID = r.attrib.get('sampleId')
        callPk = int(r.attrib.get('callPk'))
        mass = float(r.attrib.get('mass'))
        height = float(r.attrib.get('height'))
        area = float(r.attrib.get('area'))


        entry = (assayID, sampleID, callPk, mass, height, area)
#        print entry

        yield( entry )

def insert_name_string_into_assayID(assayID, name_string):
    # insert name_str before # (or at end if # is absent)
    hash_pos = assayID.find('#')
    if hash_pos == -1:
        assayID = assayID + name_string
    else:
        assayID = assayID[0:hash_pos] + name_string + assayID[hash_pos:]
    return assayID

def get_parents_and_progeny_from_file(filename, name_string):
    entries_to_exclude = set([ 'Prog_' + str(num) for num in [1, 2, 19, 32, 62, 74,75,82, 60, 126, 149, 138, 178, 216, 57, 143] ])
    assayID_sampleID_mass_height_area = get_assayID_sampleID_mass_height_area_from_XML(filename)
    for entry in assayID_sampleID_mass_height_area:
        assayID, sampleID, callPk, mass, height, area = entry

        assayID = insert_name_string_into_assayID(assayID, name_string)

        # make a new tuple with the updated assayID
        new_entry = assayID, sampleID, callPk, mass, height, area

        if sampleID == 'IACSP95-3018' or sampleID == 'IACSP93-3046':
#            print new_entry
            yield(new_entry)
        elif sampleID.find('Prog_') == 0:
            if sampleID not in entries_to_exclude:
#                print new_entry
                yield(new_entry)

def get_parents_and_progeny_from_all_files(filelist, name_string = ''):
    all_entries = []
    for filename in filelist:
        all_entries.extend( list(get_parents_and_progeny_from_file(filename, name_string)) )
    return all_entries

def build_map_from_all_files(non_batch_filelist, batch1_filelist, batch2_filelist):
    non_batch_entries = get_parents_and_progeny_from_all_files(non_batch_filelist)
    batch1_entries = get_parents_and_progeny_from_all_files(batch1_filelist, name_string = 'b1')
    batch2_entries = get_parents_and_progeny_from_all_files(batch2_filelist, name_string = 'b2')

    all_entries = non_batch_entries
    all_entries.extend(batch1_entries)
    all_entries.extend(batch2_entries)

    assayID_to_sampleID_to_data = {}
    for entry in all_entries:
        assayID, sampleID, callPk, mass, height, area = entry

        if assayID in assayID_to_sampleID_to_data:
            map_given_assayID = assayID_to_sampleID_to_data[assayID]
            if sampleID in map_given_assayID:
                if callPk in map_given_assayID[sampleID]:
                    map_given_assayID[sampleID][callPk].append( (mass, height, area) )
                else:
                    map_given_assayID[sampleID][callPk] = [(mass, height, area)]
            else:
                map_given_assayID[sampleID] = { callPk : [(mass, height, area)] }
        else:
            assayID_to_sampleID_to_data[assayID] = { sampleID : { callPk: [(mass, height, area)]} }

    for assayID in assayID_to_sampleID_to_data:
        folder = '../data/processed_data/'
        f_par = open(folder + assayID + '_parents', 'w')
        f_prog = open(folder + assayID + '_progeny', 'w')
        print 'processing locus', assayID

        map_given_assayID = assayID_to_sampleID_to_data[assayID]
        for sampleID in map_given_assayID:
            for callPk in map_given_assayID[sampleID]:
                # make sure the lower mass is the first column:
                data_point = sorted(map_given_assayID[sampleID][callPk])

                if len(data_point) != 2:
                    print 'Error: assayID', assayID, 'sampleID', sampleID, 'has', len(data_point), 'peak heights listed for callPk', callPk
                    sys.exit(1)

                low_mass, low_mass_height, low_mass_area = data_point[0]
                high_mass, high_mass_height, high_mass_area = data_point[1]
                out_string = str(sampleID) + ' ' + str(low_mass_height) + ' ' + str(high_mass_height) + ' ' + str(low_mass_area) + ' ' + str(high_mass_area) + '\n'

                if sampleID.find('Prog_') == 0:
                    f_prog.write( out_string )
#                print 'f_prog', out_string
                else:
                    f_par.write( out_string )
#                    print '\tf_par', out_string[:-1]

def find(lst, item):
	for ind,i in enumerate(lst):
		if i == item:
			return ind
	return -1

def main(argv):
    b1_index = find(argv, '--b1')
    b2_index = find(argv, '--b2')

    if b1_index == -1 or b2_index == -1 or b1_index > b2_index:
        print 'please include --b1 and --b2 in that order'
        sys.exit(1)

    non_batch_filelist = argv[:b1_index]
    batch1_filelist = argv[b1_index+1:b2_index]
    batch2_filelist = argv[b2_index+1:]

#    print 'non batch', non_batch_filelist
#    print 'batch 1', batch1_filelist
#    print 'batch 2', batch2_filelist
    build_map_from_all_files(non_batch_filelist, batch1_filelist, batch2_filelist)

if __name__ == '__main__':
    main(sys.argv[1:])
