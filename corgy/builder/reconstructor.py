def get_stem_rotation_matrix(stem, (u, v, t)):
    twist1 = stem.twists[0]

    # rotate around the stem axis to adjust the twist


    # rotate down from the twist axis
    comp1 = cross(stem.vec(), twist1)

    rot_mat1 = rotation_matrix(stem.vec(), t)
    rot_mat2 = rotation_matrix(twist1, u - pi/2)
    rot_mat3 = rotation_matrix(comp1, v)

    rot_mat4 = dot(rot_mat3, dot(rot_mat2, rot_mat1))

    return rot_mat4

def rotate_stem(stem, (u, v, t)):
    '''
    Rotate a particular stem.
    '''
    stem2 = deepcopy(stem)
    rot_mat4 = get_stem_rotation_matrix(stem, (u,v,t))
    stem2.rotate(rot_mat4, offset=stem.mids[0])

    return stem2

def rotate_chain(chain, rot_mat, offset):
    '''
    Move according to rot_mat for the position of offset.

    @param chain: A Bio.PDB.Chain instance.
    @param rot_mat: A left_multiplying rotation_matrix.
    @param offset: The position from which to do the rotation.
    '''

    atoms = Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(identity_matrix, -offset)
        atom.coord -= offset
        atom.transform(rot_mat, offset)

def translate_chain(chain, translation):
    '''
    Translate all of the atoms in a chain by a certain amount.

    @param chain: A Bio.PDB.Chain instance to be translated.
    @translation: A vector indicating the direction of the translation.
    '''
    atoms = Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(identity_matrix, translation)


def align_chain_to_stem(chain, define, stem2):
    stem1 = define_to_stem_model(chain, define)

    (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
    rot_mat = get_stem_rotation_matrix(stem1, (pi-u, -v, -t))
    rotate_chain(chain, inv(rot_mat), stem1.mids[0])
    translate_chain(chain, stem2.mids[0] - stem1.mids[0])

def reconstruct_stems(sm):
    '''
    Reconstruct the stems around a Spatial Model.

    @param sm: Spatial Model
    '''
    #sm.traverse_and_build()
    new_chain = Chain(' ')

    for stem_name in sm.stem_defs.keys():
        stem_def = sm.stem_defs[stem_name]
        stem = sm.stems[stem_name]

        filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
        #print "stem_name:", stem_name, "stem_def:", stem_def
        pdb_file = os.path.join(Configuration.stem_fragment_dir, filename)

        chain = list(PDBParser().get_structure('temp', pdb_file).get_chains())[0]
        align_chain_to_stem(chain, stem_def.define, stem)

        #print "stem_def.define:", stem_def.define
        #print "sm.bg.defines[stem_name]:", sm.bg.defines[stem_name]

        #for e in chain.get_list():
        for i in range(stem_def.bp_length+1):
            #print "i:", i
            e = chain[stem_def.define[0] + i]
            e.id = (e.id[0], sm.bg.defines[stem_name][0] + i, e.id[2])
            #print "adding:", e.id
            new_chain.add(e)

            e = chain[stem_def.define[2] + i]
            e.id = (e.id[0], sm.bg.defines[stem_name][2] + i, e.id[2])
            #print "adding:", e.id
            new_chain.add(e)

    return new_chain


def output_chain(chain, filename):
    '''
    Dump a chain to an output file.

    @param chain: The Bio.PDB.Chain to dump.
    @param filename: The place to dump it.
    '''
    m = Model(' ')
    s = Structure(' ')

    m.add(chain)
    s.add(m)

    io = PDBIO()
    io.set_structure(s)
    io.save(filename)

def define_to_stem_model(chain, define):
    '''
    Extract a StemModel from a Bio.PDB.Chain structure.

    The define is 4-tuple containing the start and end coordinates
    of the stem on each strand. 

    s1s s1e s2s s2e

    @param chain: The Bio.PDB.Chain representation of the chain
    @param define: The BulgeGraph define
    @return: A StemModel with the coordinates and orientation of the stem.
    '''
    stem = StemModel()

    mids = get_mids(chain, define)

    stem.mids = tuple([m.get_array() for m in mids])
    stem.twists = get_twists(chain, define)

    return stem

def splice_stem(chain, define):
    '''
    Extract just the defined stem from the chain and return it as
    a new chain.
    
    @param chain: A Bio.PDB.Chain containing the stem in define
    @param define: The BulgeGraph stem define
    '''
    start1 = define[0]
    end1 = define[1]

    start2 = define[2]
    end2 = define[3]

    new_chain = Chain(' ')

    for i in xrange(start1, end1+1):
        #new_chain.insert(i, chain[i])
        new_chain.add(chain[i])

    for i in xrange(start2, end2+1):
        new_chain.add(chain[i])

    '''
    m = Model(' ')
    s = Structure(' ')
    m.add(new_chain)
    s.add(m)

    io=PDBIO()
    io.set_structure(s)
    io.save('temp.pdb')
    '''

    return new_chain
