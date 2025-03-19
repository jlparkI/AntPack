/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

// C++ headers
#include <string>
#include <unordered_map>

// Library headers
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>

// Project headers
#include "annotator_classes/single_chain_annotator.h"
#include "annotator_classes/paired_chain_annotator.h"
#include "annotator_classes/annotator_base_class.h"
#include "annotator_classes/ig_aligner.h"
#include "annotator_classes/prefiltering_tool.h"
#include "vj_assignment/vj_match_counter.h"
#include "humanness_calcs/responsibility_calcs.h"
#include "dna_read_handling/dna_sequence_handling.h"
#include "utilities/utilities.h"
#include "developability_calcs/liability_search_tool.h"

namespace nb = nanobind;


NB_MODULE(antpack_cpp_ext, m) {
    nb::class_<SequenceAnnotators::AnnotatorBaseClassCpp>(m,
            "AnnotatorBaseClassCpp")
        .def(nb::init<std::string>())
        .def("sort_position_codes",
                &SequenceAnnotators::AnnotatorBaseClassCpp::sort_position_codes,
                nb::arg("position_code_list"),
     R"(
        Takes an input list of position codes for a specified scheme and
        sorts them. This is useful since for some schemes (e.g. IMGT)
        sorting is nontrivial, e.g. 112A goes before 112.

        Args:
            position_code_list (list): A list of position codes. If '-'
                is present it is filtered out and is not included in the
                returned list.

        Returns:
            sorted_codes (list): A list of sorted position codes.)")
        .def("build_msa", &SequenceAnnotators::AnnotatorBaseClassCpp::build_msa,
                nb::arg("sequences"), nb::arg("annotations"),
                nb::arg("add_unobserved_positions") = false,
     R"(
        Builds a multiple sequence alignment using a list of sequences
        and a corresponding list of tuples output by analyze_seq or
        analyze_seqs (e.g. from PairedChainAnnotator or SingleChainAnnotator).

        Args:
            sequences (list): A list of sequences.
            annotations (list): A list of tuples, each containing (numbering,
                percent_identity, chain_name, error_message). These tuples
                are what you will get as output if you pass sequences to
                the analyze_seq or analyze_seqs methods of SingleChainAnnotator
                or PairedChainAnnotator.
            add_unobserved_positions (bool): If False, only positions observed
                for one or more sequences appear in the output. If True, by
                contrast, not just observed positions but all expected positions
                for a given numbering scheme appear in the output. If IMGT expected
                position 9 does not occur in any dataset sequence, for example,
                it will not appear in the output if this argument is False but
                will be added (and will be blank for all sequences) if True.

        Returns:
            position_codes (list): A list of position codes from the appropriate numbering
                scheme.
            aligned_seqs (list): A list of strings -- the input sequences all aligned
                to form an MSA.)")
        .def("assign_cdr_labels", &SequenceAnnotators::AnnotatorBaseClassCpp::assign_cdr_labels,
                nb::arg("numbering"), nb::arg("chain"),
                nb::arg("scheme") = "",
     R"(
        Assigns a list of labels "-", "fmwk1", "cdr1", "fmwk2", "cdr2",
        "fmwk3", "cdr3", "fmwk4" to each amino acid in a sequence already
        annotated using the "analyze_seq" or "analyze_seqs" commands. The
        labels indicate which framework region or CDR each amino acid / position
        is in. This function can be used to assign CDRs with a different
        scheme than the one used to number the sequence if desired.

        Args:
            numbering (list): A list containing valid codes for the scheme that was
                selected when this object was created. If you pass a sequence to
                the analyze_seq method of SingleChainAnnotator or PairedChainAnnotator,
                the numbering will be the first element of the tuple that is returned
                (or the first element of both tuples that are returned
                for PairedChainAnnotator).
            chain (str): A valid chain (e.g. 'H', 'K', 'L', 'A'). The assigned chain is the
                third element of the tuple returned by analyze_seq. For this function
                only, 'K' and 'L' are equivalent since they both refer to a light chain,
                so if your chain is light you can supply either for the same result.
            scheme (str): Either "" or a valid scheme. If "" (default), the scheme that
                is used is the same as the one selected when the annotator was constructed.
                Using a different scheme can enable you to "cross-assign" CDRs and number
                with one scheme while assigning CDRs with another. So if you create an
                annotator with "imgt" as the scheme then you are numbering using "imgt",
                but by passing e.g. "kabat" to this function, you can use the kabat CDR
                definitions instead of the IMGT ones. Valid schemes for this
                function only are 'imgt', 'aho', 'kabat', 'martin', 'north'. For TCRs
                only "" and "imgt" are accepted.

        Returns:
            region_labels (list): A list of strings, each of which is one of
                "fmwk1", "fmwk2", "fmwk3", "fmwk4", "cdr1", "cdr2", "cdr3" or "-".
                This list will be of the same length as the input alignment.)")
        .def("trim_alignment", &SequenceAnnotators::AnnotatorBaseClassCpp::trim_alignment,
                nb::arg("sequence"), nb::arg("alignment"),
     R"(
        Takes as input a sequence and a tuple produced by
        analyze_seq and trims off any gap regions at the end
        that result when there are amino acids on either end
        which are not part of the numbered variable region.
        The output from analyze_seq can be fed directly to this
        function.

        Args:
            sequence (str): The input sequence.
            alignment (tuple): A tuple containing (numbering,
                percent_identity, chain_name, error_message). This tuple
                is what you will get as output if you pass sequences to
                the analyze_seq method of SingleChainAnnotator
                or PairedChainAnnotator.

        Returns:
            trimmed_seq (str): The trimmed input sequence.
            trimmed_numbering (list): The first element of the input tuple, the
                numbering, but with all gap regions trimmed off the end.
            exstart (int): The first untrimmed position in the input sequence.
            exend (int): The last untrimmed position in the input sequence.
                The trimmed sequence is sequence[exstart:exend].)");

    nb::class_<SequenceAnnotators::SingleChainAnnotatorCpp,
        SequenceAnnotators::AnnotatorBaseClassCpp>(m, "SingleChainAnnotatorCpp")
        .def(nb::init<std::vector<std::string>,
                std::string, std::string,
                std::unordered_map<std::string, size_t>>())
        .def("analyze_seq", &SequenceAnnotators::SingleChainAnnotatorCpp::analyze_seq,
                nb::arg("sequence"),
     R"(
        Numbers a single input sequence. A list of
        outputs from this function can be passed to build_msa
        if desired. The output from this function can also be passed
        to trim_alignment, to assign_cdr_labels and to the VJGeneTool
        as well.

        Args:
            sequence (str): A string which is a sequence containing the usual
                20 amino acids. X is also allowed but should be used sparingly.

        Returns:
            sequence_results (tuple): A tuple of (sequence_numbering, percent_identity,
                chain_name, error_message). If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.)")
        .def("analyze_seqs", &SequenceAnnotators::SingleChainAnnotatorCpp::analyze_seqs,
                nb::arg("sequences"),
     R"(
        Numbers a list of input sequences. The outputs
        can be passed to other functions like build_msa, trim_alignment,
        assign_cdr_labels and the VJGeneTool if desired.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids. X is also allowed
                (although X should be used with caution; it would be
                impossible to correctly align a sequence consisting mostly
                of X for example).

        Returns:
            sequence_results (list): A list of tuples of (sequence_numbering, percent_identity,
                chain_name, error_message). If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.)");




    nb::class_<SequenceAnnotators::PairedChainAnnotatorCpp,
        SequenceAnnotators::AnnotatorBaseClassCpp>(m, "PairedChainAnnotatorCpp")
        .def(nb::init<std::string, std::string,
                std::unordered_map<std::string, size_t>,
                std::string>())
        .def("analyze_seq", &SequenceAnnotators::PairedChainAnnotatorCpp::analyze_seq,
                nb::arg("sequence"),
     R"(
        Extracts and numbers the variable chain regions from a sequence that is
        may contain both a light ('K', 'L' for antibodies, 'B' or 'D' for TCRs)
        region and a heavy ('H' for antibodies, 'A' or 'G' for TCRs) region.
        The extracted light or heavy chains that are returned can be passed to
        other tools like build_msa, trim_alignment, assign_cdr_labels and the
        VJGeneTool.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids. X is also
                allowed but should be used sparingly.

        Returns:
            heavy_chain_result (tuple): A tuple of (numbering, percent_identity,
                chain_name, error_message). Numbering is the same length as the
                input sequence. A low percent identity or an error message may
                indicate a problem with the input sequence. The error_message is
                "" unless some error occurred.
            light_chain_result (tuple): A tuple of (numbering, percent_identity,
                chain_name, error_message). Numbering is the same length as the input
                sequence. A low percent identity or an error message may indicate a problem
                with the input sequence. The error_message is "" unless some error occurred.)")
        .def("analyze_seqs", &SequenceAnnotators::PairedChainAnnotatorCpp::analyze_seqs,
                nb::arg("sequences"),
     R"(
        Extracts and numbers the variable chain regions from a list of sequences
        may contain both a light ('K', 'L' for antibodies or 'B', 'D' for TCRs)
        region and a heavy ('H' for antibodies or 'B', 'D' for TCRs) region.
        The extracted light or heavy chains that are returned can be passed to
        other tools like build_msa, trim_alignment, assign_cdr_labels and the
        VJGeneTool.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids. X is also
                allowed but should be used sparingly.

        Returns:
            heavy_chain_results (list): A list of tuples of (numbering, percent_identity,
                chain_name, error_message). Numbering is the same length as the
                corresponding sequence. A low percent identity or an error message may
                indicate a problem with an input sequence. Each error_message is ""
                unless some error occurred for that sequence.
            light_chain_results (list): A list tuples of (numbering, percent_identity,
                chain_name, error_message). Numbering is the same length as the corresponding
                sequence. A low percent identity or an error message may indicate a problem
                with an input sequence. Each error message is "" unless some error
                occurred for that sequence.)");


    nb::class_<VJAssignment::VJMatchCounter>(m, "VJMatchCounter")
        .def(nb::init<std::map<std::string, std::vector<std::string>>,
                std::map<std::string, std::vector<std::string>>,
                nb::ndarray<double, nb::shape<22,22>, nb::device::cpu, nb::c_contig>,
                std::string, std::string>() )
        .def("assign_vj_genes",
                &VJAssignment::VJMatchCounter::assign_vj_genes,
                nb::arg("alignment"), nb::arg("sequence"),
                nb::arg("species"), nb::arg("mode") = "identity",
     R"(
        Assigns V and J genes for a sequence which has already been
        numbered, preferably by AntPack but potentially by some other
        tool. The database and numbering scheme specified when
        creating the object are used. CAUTION: Make sure the scheme
        used for numbering is the same used for the VJGeneTool.

        Args:
            alignment (tuple): A tuple containing (numbering,
                percent_identity, chain_name, error_message). This tuple
                is what you will get as output if you pass sequences to
                the analyze_seq method of SingleChainAnnotator
                or PairedChainAnnotator.
            sequence (str): A sequence containing the usual 20 amino acids -- no gaps.
                X is also allowed but should be used sparingly.
            species (str): Currently must be one of 'human', 'mouse',
                'alpaca', 'rabbit' or 'unknown'. For TCRs only 'human', 'mouse',
                'unknown' are allowed. If 'unknown', all species are checked
                to find the closest match. Note that 'unknown' will be slightly
                slower for this reason.
            mode (str): One of 'identity', 'evalue'. If 'identity' the highest
                percent identity sequence(s) are identified. If 'evalue' the
                lowest e-value (effectively best BLOSUM score) sequence(s)
                are identified.

        Returns:
            v_gene (str): The closest V-gene name(s).
            j_gene (str): The closest J-gene name(s).
            v_pident (float): If mode is 'identity', the number of positions at which
                the numbered sequence matches the v-gene divided by the total number of
                non-blank positions in the v-gene. If mode is 'evalue', the best BLOSUM
                score (this can be converted to an e-value). If more than one v-gene
                with the same score is found, multiple v-genes are returned as a single
                string delimited with '_' to separate the different v-genes.
            j_pident (float): If mode is 'identity', the number of positions at which
                the numbered sequence matches the j-gene divided by the total number of
                non-blank positions in the j-gene. If mode is 'evalue', the best BLOSUM
                score (this can be converted to an e-value). If more than one j-gene
                with the same score is found, multiple j-genes are returned as a single
                string delimited with '_' to separate the different j-genes.
            species (str): The species. This will be the same as the input species
                UNLESS your specified input species is unknown, in which case the
                species that was identified will be returned.)")
        .def("get_vj_gene_sequence",
                &VJAssignment::VJMatchCounter::get_vj_gene_sequence,
                nb::arg("vj_gene_name"), nb::arg("species"),
     R"(
        Retrieves the amino acid sequence of a specified V or J
        gene, if it is in the latest version of the specified
        database in this version of AntPack. You can use
        assign_vj_genes and this function to see what the VJ
        sequences are (if needed).

        Args:
            query_name (str): A valid V or J gene name, as generated
                by for example assign_sequence.
            species (str): One of 'human', 'mouse', 'alpaca', 'rabbit'.

        Returns:
            sequence (str): The amino acid sequence of the V or J gene
                that was requested, gapped to be length 128 consistent
                with the IMGT numbering scheme. If that V or J gene name
                does not match anything, None is returned.)")
        .def("get_seq_lists", &VJAssignment::VJMatchCounter::get_seq_lists);



    nb::class_<LiabilitySearch::LiabilitySearchToolCpp>(m,
            "LiabilitySearchTool")
        .def(nb::init<>() )
        .def("analyze_seq",
                &LiabilitySearch::LiabilitySearchToolCpp::analyze_seq,
                nb::arg("sequence"), nb::arg("alignment"),
                nb::arg("scheme"), nb::arg("cdr_scheme"),
     R"(
        Searches for some common motifs which may correspond to possible
        development liabilities. Note that this may sometimes be a false positive;
        the presence of a possible N-glycosylation motif, for example, does
        not guarantee that N-glycosylation will occur. It does however identify
        sites where there is a risk. Currently only antibodies are allowed; TCRs
        are not supported.

        Args:
            sequence (str): A sequence containing the usual 20 amino acids -- no gaps.
                X is also allowed but should be used sparingly.
            alignment (tuple): A tuple containing (numbering,
                percent_identity, chain_name, error_message). This tuple
                is what you will get as output if you pass sequences to
                the analyze_seq method of SingleChainAnnotator
                or PairedChainAnnotator.
            scheme (str): The numbering scheme. One of 'aho', 'imgt', 'kabat'
                or 'martin'. It is very important to use the same scheme
                that was used to number the sequence; using some other scheme
                may lead to incorrect motif identification.
            cdr_scheme (str): The scheme that is used for CDR definitions. This
                can be one of 'imgt', 'aho', 'martin', 'kabat', 'north'. Note
                that you can use a different set of CDR definitions than the
                numbering scheme (e.g. number with IMGT and define CDRs using
                Kabat) although usually this will be the same as 'scheme'.

        Returns:
            liabilities (list): A list of tuples. The first element of each
                tuple is a 2-tuple of (starting position, ending position)
                numbered with the start of the sequence as 0, indicating the
                start and end of the liability. The second element of each tuple
                is a string describing the type of liability found. If the list
                is empty, no liabilities were found. If sequence numbering fails,
                the list contains a single tuple indicating the cause of the failure.
                (This can occur if the expected cysteines in the chain are not present
                at the expected positions, which is also a liability.) )");



    nb::class_<DNASequenceTools::DNASeqTranslatorCpp>(m, "DNASeqTranslatorCpp")
        .def(nb::init<std::unordered_set<std::string> >() )
        .def("translate_dna_unknown_rf",
            &DNASequenceTools::DNASeqTranslatorCpp::translate_dna_unknown_rf,
            nb::arg("sequence"), nb::arg("check_reverse_complement") = false,
     R"(
        Converts a DNA sequence to protein when the reading frame is unknown
        and/or the antibody sequence may be in the forward or reverse complement.
        This function checks each possible reading frame (and if indicated the
        possible reading frames in the reverse complement) to see which is most
        likely to contain the mAb sequence based on the presence / absence
        of kmers common in mAbs and TCRs. The DNA sequence must consist
        of only A, C, T, G and N (any codon containing N will be translated
        to X, which is an allowed letter in the AntPack numbering tools) and
        should be uppercase letters. If you know which reading frame and/or
        complement the mAb sequence is in, it is generally faster to use
        translate_dna_known_rf instead, since it does not check multiple
        reading frames as this function does.

        CAUTION: The amino acid sequence generated by this function may contain
        stop codons, which are currently rejected by AntPack numbering tools.
        You may want to check output sequences for the presence of stop codons.

        Args:
            sequence (str): An uppercase sequence containing A, C, G and T.
                N is also allowed but should be used sparingly.
            check_reverse_complement (bool): If True, the reverse complement
                is also checked for a possible heavy / light chain.

        Returns:
            translated_seq (str): The function checks the possible reading frames
                and reverse complement (if indicated) for kmers common in heavy /
                light chains and returns the translated AA sequence corresponding
                to the best match found.)")
        .def("translate_dna_known_rf",
            &DNASequenceTools::DNASeqTranslatorCpp::translate_dna_known_rf,
            nb::arg("sequence"), nb::arg("reading_frame") = 0,
            nb::arg("reverse_complement") = false,
     R"(
        Converts a DNA sequence to protein when the reading frame is known
        and it is known whether the sequence is in the forward or reverse
        complement. This function translates the input sequence using the
        reading frame and complement you specify. This is faster than
        unknown_rf, but unknown_rf may be more useful in situations (e.g.
        PacBio reads) where you do not know the reading frame and forward or
        reverse complement position.

        CAUTION: The amino acid sequence generated by this function may contain
        stop codons, which are currently rejected by AntPack numbering tools.
        You may want to check output sequences for the presence of stop codons.

        Args:
            sequence (str): An uppercase sequence containing A, C, G and T.
                N is also allowed but should be used sparingly.
            reading_frame (int): One of 0, 1, or 2. Indicates how many positions
                forward to "slide" before starting translation.
            reverse_complement (bool): If True, the reverse complement is
                translated instead of the forward complement.

        Returns:
            translated_seq (str): The translated AA sequence from the input.)");


    m.def("getProbsCExt", &HumannessCalculations::getProbsCExt,
            nb::call_guard<nb::gil_scoped_release>());
    m.def("mask_terminal_deletions",
            &HumannessCalculations::mask_terminal_deletions,
            nb::call_guard<nb::gil_scoped_release>());
    m.def("getProbsCExt_masked", &HumannessCalculations::getProbsCExt_masked,
            nb::call_guard<nb::gil_scoped_release>());

    // This function is exposed for testing purposes only.
    nb::class_<PrefilteringRoutines::PrefilteringTool>(m, "PrefilteringTool")
        .def(nb::init<std::string, std::unordered_map<std::string, size_t>>())
        .def("pyfind_c_terminals",
                &PrefilteringRoutines::PrefilteringTool::pyfind_c_terminals);
}
