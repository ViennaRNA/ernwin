
import corgy.graph.bulge_graph as cgb

def define_to_group(bg, d):
    if d[0] == 's':
        return 0
    if len(bg.edges[d]) == 1:
        return 1
    if bg.weights[d] == 1:
        return 2
    return 3


def convert_graph_to_json(bg):
    '''
    Convert a BulgeGraph to JSON format to be used with the d3.js 
    force-directed graph visualization.

    @param bg: A BulgeGraph data structure
    @return: A string containing a JSON representation of the graph
    '''
    output = '{"nodes": [\n'
    bg_indeces = dict()

    counter = 0
    node_strs = []
    for d in bg.defines.keys():
        node_strs += ['{"name":"%s","group":%d}\n' % (d, define_to_group(bg, d))]
        bg_indeces[d] = counter
        counter += 1

    output += ",".join(node_strs)

    output += "],\n"
    output += '"links": ['

    link_strs = []
    visited = set()
    for d in bg.defines.keys():
        for edge in bg.edges[d]:
            if (edge,d) in visited:
                continue
            visited.add((d,edge))
            link_strs += ['{"source":%d, "target":%d, "value": 1}\n' % (bg_indeces[d],
                    bg_indeces[edge])]

    output += ",".join(link_strs)
    output += "]\n}"

    return output

def convert_graph_to_fancy_json(bg):
    '''
    Convert a BulgeGraph to JSON format to be used with the d3.js 
    force-directed graph visualization.

    The difference between this and the other graph is that in this case, the
    nodes will be the ends of the stems. The stems themselves will be edges
    as well as the bulges.

    Loop regions will be nodes connected to the end of a stem by an edge.

    @param bg: A BulgeGraph data structure
    @return: A string containing a JSON representation of the graph.
    '''
    bg_indeces = dict()

    counter = 0
    node_strs = []
    for d in bg.defines.keys():
        if d[0] == 's':
            # create names for the two ends of each stem
            d1 = "%ss0" % (d)
            d2 = "%ss1" % (d)

            node_strs += ['{"name":"%s","group":%d,"id":%d}\n' % (d1, define_to_group(bg, d), counter)]
            node_strs += ['{"name":"%s","group":%d, "id":%d}\n' % (d2, define_to_group(bg, d), counter+1)]

            bg_indeces[d1] = counter
            bg_indeces[d2] = counter+1
            counter += 2
        else:
            if len(bg.edges[d]) == 1:
                node_strs += ['{"name":"%s","group":%d,"id":%d}\n' % (d, define_to_group(bg, d), counter)]
                bg_indeces[d] = counter
                counter += 1

    link_strs = []
    distance_strs = []
    long_range_strs = []

    for d in bg.defines.keys():
        if d[0] == 's':
            d1 = "%ss0" % (d)
            d2 = "%ss1" % (d)

            link_strs += ['{"source":%d, "target":%d, "length":%d, "value": 9}\n' % (bg_indeces[d1], bg_indeces[d2], bg.defines[d][1] - bg.defines[d][0] )]
            distance_strs += ['%d' % (bg.defines[d][1] - bg.defines[d][0])]
            continue
        
        if len(bg.edges[d]) == 1:
            d_bulge = d
            d_stem = list(bg.edges[d])[0]

            (s1b, s1e) = bg.get_sides(d_stem, d_bulge)

            d1 = "%ss%d" % (d_stem, s1b)
            d2 = d_bulge

            link_strs += ['{"source":%d, "target":%d, "length":%d, "value": 1}\n' % (bg_indeces[d1], bg_indeces[d2], bg.defines[d][1] - bg.defines[d][0])]
            distance_strs += ['%d' % (bg.defines[d][1] - bg.defines[d][0])]
        else:
            bgl = list(bg.edges[d])
            s1 = bgl[0]
            s2 = bgl[1]

            (s1b, s1e) = bg.get_sides(s1, d)
            (s2b, s2e) = bg.get_sides(s2, d)

            d1 = "%ss%d" % (s1, s1b)
            d2 = "%ss%d" % (s2, s2b)

            l1 = bg.defines[d][1] - bg.defines[d][0]
            if bg.weights[d] == 2:
                l2 = bg.defines[d][3] - bg.defines[d][2]
            else:
                l2 = l1

            link_strs += ['{"source":%d, "target":%d, "length":%d, "value": 1}\n' % (bg_indeces[d1], bg_indeces[d2], min(l1,l2))]
            distance_strs += ['%d' % (min(l1, l2))]

    visited = set()
    for key1 in bg.longrange.keys():
        for key2 in bg.longrange[key1]:
            if (key2, key1) in visited:
                continue
            visited.add((key1, key2))
            ends = []
            for key in [key1, key2]:
                # if connection a stem, then the connection will
                # be halfway between the start and end
                if key[0] == 's':
                    ends += [bg_indeces[key + 's0'], bg_indeces[key + 's1']]
                else:
                    if len(bg.edges[key]) == 1:
                        ends += [bg_indeces[key], bg_indeces[key]]
                    else:
                        bgl = list(bg.edges[key])

                        s1 = bgl[0]
                        s2 = bgl[1]

                        (s1b, s1e) = bg.get_sides(s1, key)
                        (s2b, s2e) = bg.get_sides(s2, key)

                        d1 = "%ss%d" % (s1, s1b)
                        d2 = "%ss%d" % (s2, s2b)

                        ends += [bg_indeces[d1], bg_indeces[d2]]

            
            long_range_strs += ['{"from1": %d, "from2":%d, "to1": %d, "to2":%d}\n' % (ends[0], ends[1], ends[2], ends[3])]

    output = '{"nodes": [\n'
    output += ",".join(node_strs)
    output += "],\n"
    '''
    output += '"distances": ['
    output += ",".join(distance_strs)
    output += "],\n"
    '''
    output += '"long_range": ['
    output += ",".join(long_range_strs)
    output += "],\n"

    output += '"links": ['
    output += ",".join(link_strs)
    output += "]\n}"

    return output
