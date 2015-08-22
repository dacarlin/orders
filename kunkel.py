import math
import json
import csv
import datetime
from autoprotocol.util import make_dottable_dict
from autoprotocol.pipette_tools import *
from autoprotocol.container import Container, WellGroup
from autoprotocol.protocol import Ref
from autoprotocol.unit import Unit


import transcriptic

config = transcriptic.Config.from_file('~/.transcriptic')


def return_media(mtype):
    '''
        Dicts of all plates available that can be purchased.
    '''
    if mtype == 'solid':
        media = {"50_ugml_Kanamycin": "ki17rs7j799zc2",
                 "100_ugml_Ampicillin": "ki17sbb845ssx9",
                 "100_ugmL_Spectinomycin": "ki17sbb9r7jf98",
                 "noAB": "ki17reefwqq3sq"}
    elif mtype == 'liquid':
        media = {"50_ugml_Kanamycin": "lb-broth-50ug-ml-kan",
                 "100_ugml_Ampicillin": "lb-broth-100ug-ml-amp",
                 "100_ugmL_Spectinomycin": "lb-broth-100ug-ml-specto"}

    else:
        raise ValueError("mtype has to be solid or liquid")
    return(media)


def det_new_group(i, base=0):
    '''
        Helper to determine if new_group should be added. Returns true when i matches the base, which defaults to 0.
    '''
    assert isinstance(i, int), "Needs an integer."
    assert isinstance(base, int), "Base has to be an integer"
    if i == base:
        new_group = True
    else:
        new_group = False
    return new_group


def printdatetime(time=True):
    printdate = datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    if not time:
        printdate = datetime.datetime.now().strftime('%Y-%m-%d')
    return printdate


def provision_to_tube(protocol, name, tube, resource_id, volume, discard=True, storage=None):
    '''
        provision_to_tube allows us to provision into a tube/well that can then be used
        in transfers on the workcell.
    '''
    assert isinstance(volume, (Unit, int, float)), "Volume must be of type int, float or Unit."
    if isinstance(volume, Unit):
        volume = volume.value
    if storage:
        dest = protocol.ref(name, None, tube, storage=storage).well(0)
    else:
        dest = protocol.ref(name, None, tube, discard=discard).well(0)
    protocol.provision(resource_id, dest, "%s:microliter" % volume)
    return(dest)


def transfer_kwargs(pre_buffer, one_tip=False, one_source=False):
    kwargs = {"one_tip": one_tip,
              "one_source": one_source,
              "pre_buffer": "%s:microliter" % pre_buffer,
              "blowout_buffer": True}
    return(kwargs)


def thermocycle_ramp(start_temp, end_temp, total_duration, step_duration):
    '''
        Create a ramp instruction for the thermocyler. Used in annealing protocols.
    '''
    assert Unit.fromstring(total_duration).unit == Unit.fromstring(step_duration).unit, ("Thermocycle_ramp durations"
                                                                                         " must be specified using the"
                                                                                         " same unit of time.")
    thermocycle_steps = []
    start_temp = Unit.fromstring(start_temp).value
    num_steps = int(Unit.fromstring(total_duration).value // Unit.fromstring(step_duration).value)
    step_size = (Unit.fromstring(end_temp).value - start_temp) // num_steps
    for i in xrange(0, num_steps):
        thermocycle_steps.append({
            "temperature": "%d:celsius" % (start_temp + i * step_size), "duration": step_duration
            })
    return thermocycle_steps


def ref_kit_container(protocol, name, container, kit_id, discard=True, store=None):
    '''
        Still in use to allow booking of agar plates on the fly
    '''
    kit_item = Container(None, protocol.container_type(container))
    if store:
        protocol.refs[name] = Ref(name, {"reserve": kit_id, "store": {"where": store}}, kit_item)
    else:
        protocol.refs[name] = Ref(name, {"reserve": kit_id, "discard": discard}, kit_item)
    return(kit_item)


def find_aliquot_by_name(name):
    # looks up code for the sample
    resp = config.get('samples?q=%s&show_containers=false' % name)
    for aq in resp.json()['results']:
        if aq['name'] == name:
            return aq


def kunkel_mutagenesis(protocol, params):
    def find_part(part_name, check_primer=False):
        if check_primer and (params.t7pro or params.t7term):
            return 'primer'
        else:
            aq = find_aliquot_by_name(part_name)
            if not aq:
                raise ValueError(
                    "Couldn't find aliquot with name '%s'" % part_name)
            container = None
            for r in protocol.refs.itervalues():
                if r.opts.get('id') == aq['container']['id']:
                    container = r.container
                    break
            if container is None:
                container = protocol.ref(
                    aq['container']['label'] or aq['container']['id'],
                    aq['container']['id'],
                    aq['container']['container_type']['shortname'],
                    discard=False,
                    storage="cold_20"
                    )
            return container.well(aq['well_idx']).set_volume("%s:microliter" % aq['volume_ul'])

    def make_10_atp(vol):
        # account for tube dead_vol
        vol = vol + 15
        atp = protocol.ref("atp_10mM", cont_type='micro-1.5', discard=True).well(0)
        protocol.provision('rs16pccshb6cb4', atp, "%s:microliter" % (vol/10))
        protocol.provision('rs17gmh5wafm5p', atp, "%s:microliter" % (vol - (vol/10)))
        return atp

    def provision_reagents(reagents, dest):
        for key, reagent in reagents.iteritems():
            protocol.provision(reagent[0],
                               dest,
                               "%s:microliter" % ((len(params.mutants) + 1.0) * reagent[1] * add_mm))

    def isLast(itr):
        old = itr.next()
        for new in itr:
            yield False, old
            old = new
        yield True, old

    params = make_dottable_dict(params)
    params.ssDNA = find_part(params.ssDNA)
    params.mutants = []
    # read in oligos and mutants
    with open('kunkel_mutants.csv', 'rU') as rxtmap:
        reader = csv.reader(rxtmap)
        current_mutant_label = ''
        current_sequencing_primer = ''

        for idx, (is_last, row) in enumerate(isLast(reader)):
            if idx == 0 and row[0] == 'mutant_label':
                continue
            if row[0] != current_mutant_label:
                if current_mutant_label != '':
                    params.mutants.append({
                        'sequencing_primer': find_part(current_sequencing_primer, True),
                        'mutant_label': current_mutant_label,
                        'oligos': my_oligos
                    })
                current_mutant_label = row[0]
                current_sequencing_primer = row[1]
                my_oligos = []

            my_oligos.append({"oligo_label": row[2],
                              "sequence": row[3],
                              "scale": row[4],
                              "purification": row[5]})
            if is_last:
                params.mutants.append({
                        'sequencing_primer': find_part(current_sequencing_primer, True),
                        'mutant_label': current_mutant_label,
                        'oligos': my_oligos
                    })

    # mastermix vol to make - needs to be adjusted based on testing
    add_mm = 1.3
    # Get unique set of oligos based on sequence data
    # Oligosynthesize
    oligos_to_synthesize = []
    for i, mutant in enumerate(params.mutants):
        mutant["mutant_label"] = mutant["mutant_label"] or "mutant_%s" % (i + 1)
        oligos_to_synthesize.append(mutant["oligos"])
    flattened = [val for sublist in oligos_to_synthesize for val in sublist]
    oligos_to_synthesize = list({v['sequence']: v for v in flattened}.values())

    # re-factor to remove add_properites
    oligo_containers = []
    oligos = []
    for i, oligo in enumerate(oligos_to_synthesize):
        label = oligo["oligo_label"] or "seq_%s" % i
        oligo_containers.append(protocol.ref(label, None, "micro-2.0",
                                storage="cold_4").well(0))
        oligo_containers[i].add_properties({"sequence": oligo["sequence"]})
        oligos.append({"sequence": oligo["sequence"],
                       "destination": oligo_containers[i],
                       "scale": oligo["scale"],
                       "purification": oligo["purification"]})

    protocol.oligosynthesize(oligos)

    # Kinase
    kinase_oligo_plate = protocol.ref("kinase_oligo_plate", None, "96-pcr",
                                      storage="cold_20")
    wells_to_kinase = kinase_oligo_plate.wells_from("A1", len(oligos))

    # provision atp for entire protocol
    atp_needed = ((len(oligos) * add_mm) + (len(params.mutants) + 1) * 0.4 * add_mm)
    atp = make_10_atp(atp_needed)

    kinase_mix = []
    for i in range(int(math.ceil(len(oligos)/60.0))):
        kinase_mix.append(protocol.ref("kinase_mix-%s" % (i + 1), None, "micro-1.5", discard=True).well(0))
    reagents = {'pnkbuffer': ['rs16pc9rd5sg5d', 3],
                'water': ['rs17gmh5wafm5p', 18],
                'pnk': ['rs16pc9rd5hsf6', 1]}
    provision_reagents(reagents, kinase_mix)
    protocol.transfer(atp, kinase_mix, "%s:microliter" % ((len(params.mutants) + 1) * 1 * add_mm), new_group=True)

    protocol.transfer(kinase_mix,
                      wells_to_kinase,
                      "23:microliter",
                      **transfer_kwargs(15, True, True))

    for i, oligo in enumerate(oligo_containers):
        protocol.transfer(oligo,
                          wells_to_kinase[i],
                          "7:microliter",
                          mix_after=False,
                          new_group=det_new_group(i),
                          aspirate_source=aspirate_source(depth=depth("ll_following",
                                                          lld="pressure",
                                                          distance="0.0:meter")),
                          **transfer_kwargs(10))

    protocol.seal(kinase_oligo_plate)

    protocol.thermocycle(kinase_oligo_plate,
                         [{"cycles": 1, "steps": [
                             {"temperature": "37:celsius",
                              "duration": "60:minute"},
                             ]}
                          ], volume="30:microliter")

    # make ssDNA_mastermix
    mix_plate = protocol.ref("mix_plate", None, "96-pcr", discard=True)
    ssDNA_mix = mix_plate.well(0)
    protocol.transfer(params.ssDNA,
                      ssDNA_mix,
                      "%s:microliter" % ((len(params.mutants) + 1) * 2.0 * add_mm),
                      **transfer_kwargs((len(params.mutants) + 1) * 1))
    protocol.provision('rs17sh5rzz79ct', ssDNA_mix, "%s:microliter" % ((len(params.mutants) + 1) * 0.2 * add_mm))

    # Dilute
    protocol.unseal(kinase_oligo_plate)

    diluted_oligo_plate = protocol.ref("dilute_oligo_plate", None, "96-flat", discard=True)
    diluted_oligo_wells = diluted_oligo_plate.wells_from(0, len(params.mutants))

    water = [provision_to_tube(protocol, "water%s" % (i + 1), "micro-2.0", "rs17gmh5wafm5p", 1900)
             for i in range(int(math.ceil(len(params.mutants)/float(9.5))))
             ]

    protocol.transfer(water,
                      diluted_oligo_wells,
                      "200:microliter",
                      disposal_vol="0:microliter",
                      **transfer_kwargs(40, True, True))

    mutants = [m for m in params.mutants if m]
    mutants = sorted(mutants, key=lambda mutant: mutant["mutant_label"])

    for i, m in enumerate(mutants):
        for j, kin_oligo in enumerate(m["oligos"]):
            if i == 0 and j == 0:
                new_group = True
            else:
                new_group = False
            index = next((i for i, olig in enumerate(oligo_containers) if olig.properties["sequence"] == kin_oligo["sequence"]), -1)
            protocol.transfer(kinase_oligo_plate.well(index),
                              diluted_oligo_plate.well(i),
                              "2:microliter",
                              mix_after=True,
                              mix_vol="2:microliter",
                              new_group=new_group,
                              **transfer_kwargs(10))

    protocol.cover(diluted_oligo_plate)
    protocol.spin(diluted_oligo_plate, "700:meter/second^2", "2:minute")
    protocol.uncover(diluted_oligo_plate)

    # Anneal

    annealing_plate = protocol.ref("annealing_oligo_plate", None, "384-pcr", storage="cold_20")
    anneal_wells = annealing_plate.wells_from(0, len(params.mutants))

    protocol.transfer(ssDNA_mix,
                      anneal_wells.wells,
                      "2.2:microliter",
                      dispense_speed="50:microliter/second",
                      **transfer_kwargs(7, True, True))

    for i, oligo_reaction in enumerate(zip(diluted_oligo_wells.wells, anneal_wells.wells)):
        protocol.transfer(oligo_reaction[0],
                          oligo_reaction[1],
                          "2:microliter",
                          aspirate_source=aspirate_source(depth("ll_bottom", distance=".001:meter")),
                          mix_after=True,
                          mix_vol="2:microliter",
                          flowrate="50:microliter/second",
                          repetitions=2,
                          new_group=det_new_group(i),
                          **transfer_kwargs(5))

    protocol.seal(annealing_plate)
    protocol.spin(annealing_plate, "700:meter/second^2", "2:minute")
    protocol.thermocycle(annealing_plate, [{
        "cycles": 1,
        "steps": thermocycle_ramp("95:celsius", "25:celsius", "60:minute", "4:minute")
        }],
        volume="5:microliter",
        dataref=None,
        dyes=None)

    # Step 4 - Polymerize

    protocol.unseal(annealing_plate)
    polymerize_MM = mix_plate.well(12)
    reagents = {"buffer": ['rs17sh5rzz79ct', 0.6],
                "t4ligase": ['rs16pc8krr6ag7', 0.4],
                "t7polymerase": ['rs16pca2urcz74', 0.4],
                "dntp": ['rs16pcb542c5rd', 0.4]
                }
    provision_reagents(reagents, polymerize_MM)
    protocol.transfer(atp, polymerize_MM, "%s:microliter" % ((len(params.mutants) + 1) * 0.4 * add_mm), new_group=True)

    for reaction in anneal_wells.wells:
        protocol.transfer(polymerize_MM,
                          reaction,
                          "2.2:microliter",
                          mix_after=False,
                          **transfer_kwargs(10))

    protocol.seal(annealing_plate)
    protocol.incubate(annealing_plate, "ambient", "1.5:hour")

    # Transformation using Zymo 10B Competent Cells
    transformation_cells = []
    for i in range(len(params["mutants"])):
        transformation_cells.append(provision_to_tube(protocol, "cell_%s" % (i), "micro-1.5",
                                                      "rs16pbjc4r7vvz", 50))

    num_colonies = params["num_colonies"]
    assert len(params["mutants"]) * num_colonies <= 96, ("This protocol is limited to 96 sequenced colonies, please "
                                                         "submit additional runs if needed.")
    transformation_plate = protocol.ref("transformation_plate", None, "96-pcr", discard=True)
    protocol.incubate(transformation_plate, "cold_20", "10:minute")
    transformation_wells = transformation_plate.wells_from(0, len(params.mutants), columnwise=False)

    for i, tube in enumerate(transformation_cells):
        protocol.transfer(tube, transformation_wells[i], "50:microliter", mix_after=False)

    protocol.unseal(annealing_plate)

    for i, rxt in enumerate(anneal_wells):
        protocol.transfer(anneal_wells[i], transformation_wells[i], "2.0:microliter",
                          dispense_speed="10:microliter/second",
                          mix_after=False,
                          new_group=det_new_group(i))

    protocol.cover(transformation_plate, lid="universal")
    protocol.incubate(transformation_plate, "cold_4", "20:minute", shaking=False, co2=0)
    protocol.uncover(transformation_plate)

    agar_plates = []
    assert len(mutants) == len(transformation_wells), ("Sanity check failed. There is an issue with the number of"
                                                       "mutants and the number of transformations.")
    for well in range(0, len(transformation_wells), 6):
        agar_plate = ref_kit_container(protocol,
                                       "agar-%s_%d_%s" % (params["antibiotic"].split("_")[-1],
                                                          well + 1,
                                                          printdatetime(time=False)),
                                       "6-flat",
                                       return_media('solid')[params["antibiotic"]],
                                       discard=False, store='cold_4')
        agar_plates.append(agar_plate)
        for i, well in enumerate(transformation_wells[well:well + 6]):
            protocol.spread(well, agar_plate.well(i), "50:microliter")
        protocol.incubate(agar_plate, "warm_37", "18:hour")

    growth_plate = protocol.ref("growth_plate_%s" % printdatetime(time=False), None, "96-flat", discard=True)

    cols = int(math.ceil(len(params.mutants) * num_colonies / float(8)))
    columns = [{"column": i, "volume": "150:microliter"} for i in range(0, cols)]

    protocol.dispense(growth_plate, return_media('liquid')[params["antibiotic"]], columns)

    growth_wells = growth_plate.wells_from(0, num_colonies*len(params.mutants), columnwise=True)

    i = 0
    for k, plate in enumerate(agar_plates):
        for j in range(6):
            if plate.well(j).volume:
                protocol.autopick(plate.well(j),
                                  growth_wells[i:i+num_colonies],
                                  min_count=1,
                                  dataref=mutants[k * 6 + j]['mutant_label'])
                i = i + num_colonies

    protocol.cover(growth_plate, lid="low_evaporation")
    protocol.incubate(growth_plate, "warm_37", "24:hour", shaking=True, co2=0)
    protocol.uncover(growth_plate)

    if params.t7pro:
        seq_plate = protocol.ref("sequencing_plate_t7pro_%s" % printdatetime(time=False),
                                 cont_type="96-pcr",
                                 storage="cold_4")
        seq_well_group = seq_plate.wells_from(0, num_colonies*len(params.mutants), columnwise=True)
        t7pro_primer = protocol.ref("t7promoter",
                                    cont_type="micro-1.5",
                                    storage="cold_4")
        protocol.provision("rs17tcpekfy7v9", t7pro_primer.well(0), "1:microliter")
        protocol.provision("rs17gmh5wafm5p", t7pro_primer.well(0), "%s:microliter" % (num_colonies*(len(params.mutants)+2)))
        protocol.transfer(growth_wells, seq_well_group, "30:microliter")
        protocol.seal(seq_plate)
        protocol.sangerseq(seq_plate, seq_well_group.indices(), "Seq_primer_T7promoter", type="rca", primer=t7pro_primer)

    if params.t7term:
        seq_plate = protocol.ref("sequencing_plate_t7term_%s" % printdatetime(time=False),
                                 cont_type="96-pcr",
                                 storage="cold_4")
        seq_well_group = seq_plate.wells_from(0, num_colonies*len(params.mutants), columnwise=True)
        t7term_primer = protocol.ref("t7terminator",
                                    cont_type="micro-1.5",
                                    storage="cold_4")
        protocol.provision("rs17tcpwfbgzqd", t7term_primer.well(0), "1:microliter")
        protocol.provision("rs17gmh5wafm5p", t7term_primer.well(0), "%s:microliter" % (num_colonies*(len(params.mutants)+2)))
        protocol.transfer(growth_wells, seq_well_group, "30:microliter")
        protocol.seal(seq_plate)
        protocol.sangerseq(seq_plate, seq_well_group.indices(), "Seq_primer_t7terminator", type="rca", primer=t7term_primer)

    if not (params.t7pro or params.t7term):
        seq_plate = protocol.ref("sequencing_plate_%s" % printdatetime(time=False),
                                 cont_type="96-pcr",
                                 storage="cold_4")
        seq_well_group = seq_plate.wells_from(0, num_colonies*len(params.mutants), columnwise=True)
        j = 0
        seq_primers = {}
        mutant_well_table = {}
        for seq_set in mutants:
            seq_wells = WellGroup(seq_well_group[j:j+num_colonies])
            protocol.transfer(growth_wells[j:j+num_colonies], seq_wells, "30:microliter")
            if seq_set["sequencing_primer"] not in seq_primers:
                seq_primers[seq_set["sequencing_primer"]] = seq_wells.indices()
                mutant_well_table[seq_set["sequencing_primer"]] = {seq_set["mutant_label"]: seq_wells.indices()}
            else:
                seq_primers[seq_set["sequencing_primer"]].extend(seq_wells.indices())
                if seq_set["mutant_label"] not in mutant_well_table[seq_set["sequencing_primer"]]:
                    mutant_well_table[seq_set["sequencing_primer"]].update({seq_set["mutant_label"]: seq_wells.indices()})
                else:
                    mutant_well_table[seq_set["sequencing_primer"]][seq_set["mutant_label"]].extend(seq_wells.indices())
            j += num_colonies

        protocol.seal(seq_plate)

        for primer, wells in seq_primers.iteritems():
            dataref = "Seq_primer_%s" % (primer.container.name)
            assert primer.volume - Unit("1.1", "microliter").__mul__(len(wells)) > Unit("0", "microliter"), "You must have at least 1.1uL of sequencing primer per reaction well."
            protocol.sangerseq(seq_plate, wells, dataref, type="rca", primer=primer.container)

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(kunkel_mutagenesis)
