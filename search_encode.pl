#!/usr/bin/perl

# search_encode.pl: A general framework to query and retrieve data from the
# ENCODE repository using the REST API. Copyright (C) 2015, Adam Diehl,
# adadiehl@umich.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use JSON;
use WWW::Mechanize;
use Getopt::Long;
use Encode qw(encode_utf8 decode_utf8);
use Digest::MD5 qw(md5_hex);
use URI::Escape;
use File::Which;
#use Devel::Size qw(size total_size);

my $usage =
"\nsearch_encode.pl

A general framework to query and retrieve experiment data from the ENCODE
repository. Capable of returning metadata for search results in tabular format
or as json data, as well as downloading all or user-specified types of
datafiles from the resulting experiments.


USAGE:

search_encode.pl <search term> <params string> [OPTIONS]

<search term>
    Text string to search for within the ENCODE portal.

<params string>
    Comma-delimited sting of key:value pairs corresponding to columns in the
    given database table. See the schema at the address above. Some columns
    have sub-keys that can be searched. Some examples are available in Box 3:
    API Resources, in \"Deciphering ENCODE\", available at
    http://www.cell.com/trends/genetics/fulltext/S0168-9525%2816%2900017-2
    
    Possible keys of interest for experiments:
        biosample_term_name: tissue/cell used in assay
        assay_term_name: type of assay (e.g., ChIP-seq)
        target.investigated_as: type of target (e.g, transcription factor)
        target.label: name of factor probed (e.g., CTCF)
        status: status of experiment (released/revoked)

    Use the empty string \"\" for no parameters.

    NOTE: Due to changes in the ENCODE portal after the official switch
    to hg38 as the default human genome assembly, it is now advisable to
    specify the desired genome assembly (hg19 or hg38) for file downloads
    as a json parameter, using '--filter-json \"assembly=hg19\"', for
    example. See the section below for the --filter-json option for more
    details. 


OPTIONS:

--download
    Download files associated with the results. Type(s) of files may be
    specified with the --file-type flag. Currently only tested for results
    of type \"experiment\". May behave unpredictably for other types!

--file-list
    Instead of downloading files, print a table of metadata for files
    matching the given criteria. Compatible with the same options as
    --download.

--quality-data
    Compile a table of quality metrics for experiments matching the query,
    to include read count, % mapped reads, % duplication, NSC, RSC, cc(fragLen),
    cc(phantomPeak), cc(min), PBC, and NRF. See definitions of these metrics
    here:
    https://genome.ucsc.edu/ENCODE/qualityMetrics.html#definitions
    and in the ENCODE ChIP-seq guidelines:
    Landt, S. G., et al. (2012). ChIP-seq guidelines and practices of the
    ENCODE and modENCODE consortia. Genome Res, 22(9), 1813–1831.

--reads-and-controls
    Download sequencing reads and control reads for matching experiments.
    Forces --download and --output-type <reads>.

--output-type <type_1,...,type_n>
    Used with --download, limits results to files of given type(s).
    Multiple types may be supplied as a comma-separated list. File types are
    matched against the \"output_type\" column of the file records, and
    available values will vary depending on the experiment. Common values are:

    * reads: raw sequnence reads (usually fastq format)
    * alignments: sequence alignment (usually BAM format)
    * peaks: processed peaks from some peak caller (often bed or bigBed)
    * signal: raw signal (often bigWig format)

--file-format <format_1,...,format_n>
    Used with --output-type, restrict downloads to a given file format, e.g.,
    bam, bigBed. Available values vary depending on the specific experiment.
    This MUST contain a value for every file type supplied with --output-type
    or it will throw an error. Specify \"0\" for all types not requiring a
    file type.
    
--file-format-type <type_1,...,type_n>
    Used with --output-type, restrict downloads to a given file format type,
    e.g., \"broadPeak\" or \"narrowPeak\". See notes for --file-format.

--download-path <path>
    Location (absolute or relative to workding directory) to write downloaded
    files.

--out-root <root>
    String to prepend to output files.

--check-exists
    Check to see if a file already exists before downloading and do not
    retrieve if found.

--use-wget
    Retrieve files with the system command \"wget\". This can be faster
    for large result sets. Wget must be installed on the system for this
    to work. Checksums will not be compared with this method.

--save-json
    Save the JSON metadata for every file downloaded. This can be useful for
    identifying individual files and properties later on when retrieving
    multiple files, especially if further programmatic processing is
    anticipated.

--no-header
    Do not print a header line with column names to the metadata file.

--count-only
    Return only the number of results for the search. Do not retrieve
    any metadata or files. Suppresses --download (and related options) and
    --save-json.

--by-biosample
    Convert the search term to a biosample term name, for which search
    results will be returned. The search term is not directly queried with
    this option, but rather converted into a search parameter. This can be
    useful when a search for a specific cell line is yielding either no
    results or imprecise results.

--filter-json \"key1=val1,...,keyN=valN\"
    Search the JSON data for matching key-value pairs not available through
    the ENCODE search API. For nested terms, use format x.y.z=val, where
    x, y and z are the nested json attributes. Works for both hash elements
    (enclosed in curly brackets in the json data: {...}), and array elements
    (enclosed in square brackets in the json data: [...]), but with slightly
    different behavior. For hash elements, the value stored within the hash
    key given must be an exact match to the value. For array elements, a
    record will be accepted if ANY values within the array exactlty match
    the query value. It is currently only possible to search within named
    elements, thus hash/array elements nested within another array are not
    available.

--exclude-json \"key1=val1,...,keyN=valN\"
    Like --filter-json, except that files with any matching key-value pairs
    in their json records will be excluded from results.

--match-all
    Used with --filter-json, only return records matching all filters.

--match-all-excl
    Used with --exclude-json, only exclude records matching all filters.

--random <N>
    Select <N> results at random.

--help
    Display this message


NOTES:

* A search string of \"*\" will return all results matching the rest of the
parameters. This can be useful if you want to retrieve the same type of data
for all available cell types, for instance.

* Most table columns in the ENCODE database tables are visible to the search
algorithm and you may chain as many as you would like onto your query.
Please consult the schema, at https://www.encodeproject.org/profiles/graph.svg
to see the columns for each data table.

* Columns in tables referenced by the primary table may also be searched. (All
links are described in the schema) For instance, to search for experiments for
human only based on scientific name, you need to use the \"replicates\"
coumn, which (indirectly) references the \"organism\" table:
&replicates.library.biosample.donor.organism.scientific_name=Homo sapiens

* Fields in the JSON with values specified as a folder location usually are
available to directly search through the parameter string by specifying as
&field.column=value. E.g., &documents.attachment.download=file.pdf will
find files referencing file.pdf.

* All column names and values are CaSe-SeNsItIvE

* Database values may include white-space. To include a term with white space
in your query, put it in quotes when supplying it at the command line.

* The same key may be used multiple times with different values and results
matching any of the values will be returned. This can be useful to retrieve
data for multiple cell types or transcription factors, for example.

* An excellent tutorial on using the API and navigating the JSON data is
available at: https://www.encodeproject.org/help/rest-api/


EXAMPLES

* Find all experiments that used the GAIIx protocol from Michael Snyder's
lab and print results as a data table:

search_encode.pl \"*\" \"&documents.attachment.download=ChIP-Seq_protocol_Snyder_lab_GAIIx.pdf\"

* Same as above, but only for experiments performed on GM12878 cells:

search_encode.pl \"GM12878\" \"&documents.attachment.download=ChIP-Seq_protocol_Snyder_lab_GAIIx.pdf\"

* Same as previous, but download all peak annotations in bigBed format:

search_encode.pl \"GM12878\" \"&documents.attachment.download=ChIP-Seq_protocol_Snyder_lab_GAIIx.pdf\" --download --output-type \"peaks\" --file-format \"bigBed\"


CREDITS AND LICENSE:

Copyright (C) 2016, Adam Diehl

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Any questions/comments regarding this program may be directed to Adam Diehl:
adadiehl\@umich.edu

Originally published in \"Deciphering ENCODE\". If you use this tool in your
research, please cite:

Diehl AD and Boyle AP (2016) Deciphering ENCODE.
Trends in Genetics. Volume 32, Issue 4, 238-249.

\n";

my $search_str = $ARGV[0];
my $params = $ARGV[1];

my $help = 0;
my $download = 0;
my $file_list = 0;
my $quality_data = 0;
my $output_type_str;
my $file_format_str;
my $file_format_type_str;
my $download_path;
my $out_root;
my $save_json = 0;
my $rec_type = "Experiment";
my $no_header = 0;
my $count_only = 0;
my $by_biosample = 0;
my $check_exists = 0;
my $use_wget = 0;
my $filter_json_str;
my $excl_json_str;
my $match_all = 0;
my $match_all_excl = 0;
my $debug = 0;
my $n_rand = 0;
my $is_reads = 0;
my $reads_and_controls = 0;

GetOptions (
    "help" => \$help,
    "download" => \$download,
    "file-list" => \$file_list,
    "quality-data" => \$quality_data,
    "reads-and-controls" => \$reads_and_controls,
    "output-type=s" => \$output_type_str,
    "file-format=s" => \$file_format_str,
    "file-format-type=s" => \$file_format_type_str,
    "download-path=s" => \$download_path,
    "out-root=s" => \$out_root,
    "save-json" => \$save_json,
    "rec-type=s" => \$rec_type,
    "no-header" => \$no_header,
    "count-only" => \$count_only,
    "by-biosample" => \$by_biosample,
    "check-exists" => \$check_exists,
    "use-wget" => \$use_wget,
    "filter-json=s" => \$filter_json_str,
    "exclude-json=s" => \$excl_json_str,
    "match-all" => \$match_all,
    "match-all-excl" => \$match_all_excl,
    "random=i" => \$n_rand,
    "debug" => \$debug
    );

# Check for proper usage and help option and exit with usage message as needed.
if ($help || $#ARGV+1 < 2) {
    die "$usage\n";
}

# Check for other options and do some sanity checks...
if (defined($output_type_str) && !($download || $file_list || $quality_data) ) {
    die "\n--output-type only works with --download or --file-list. Try --help.\n";
}
#if (defined($file_format_type_str) && !defined($file_format_str)) {
#    die "\n--file-format-type requires --file-format. Try --help.\n";
#}
if (defined($download_path) && !$download) {
    die "\n--download-path requires --download. Try --help\n";
}
if ($download && !defined($output_type_str)) {
    print STDERR "\nWarning: --download will retrieve all accepted files associated with search results. If this is not what you want, use --file-type to specify which type(s) you would like to download. See --help.\n\n";
}
if ($rec_type ne "Experiment") {
    print STDERR "\nWarning: --rec_type other than \"Experiment\" has not been well tested and may produces unpredictable results!\n";
}
if ($count_only && $download) {
    print STDERR "\nWarning: --count-only overrides --download. No data will be retrieved! (see --help if this is not what you want)\n";
}
if ($use_wget) {
    my $wg_path = which('wget');
    if (!-e $wg_path) {
	die "Cannot find wget on your system. Please try again without --use-wget.\n";
    }
}
if (defined($filter_json_str) && !($download || $file_list) && !$quality_data) {
    print STDERR "\nWarning: --filter-json has no effect without --download or --file-list.\n";
}
if (defined($filter_json_str) && !($download || $file_list) && $quality_data) {
    print STDERR "\nWarning: --filter-json invoked with --quality-data. File-level filtering will be used but no files will be downloaded. Experiment-level quality metrics for experiments with files matching criteria will be written to metadata.\n";
}
if (defined($excl_json_str) && !($download || $file_list) && !$quality_data) {
    print STDERR "\nWarning: --exclude-json has no effect without --download or --file-list.\n";
}
if (defined($excl_json_str) && !($download || $file_list) && $quality_data) {
    print STDERR "\nWarning: --exclude-json invoked with --quality-data. File-level filtering will be used but no files will be downloaded. Experiment-level quality metrics for experiments with files matching criteria will be written to metadata.\n";
}

my %controls;
if ($reads_and_controls) {
    $is_reads = 1;
    $output_type_str = "reads";
    $download = 1;
}

# Check the download path for a trailing slash and add one if needed
if (defined($download_path)) {
    if ($download_path !~ m/\/$/) {
	$download_path .= '/';
    }
}

# If file types have been specified, process the string(s).
my @output_types;
if (defined($output_type_str)) {
    @output_types = split /,/, $output_type_str;
    foreach my $type (@output_types) {
	if ($type eq "reads") {
	    $is_reads = 1;
	}
    }
}


my @file_formats;
if (defined($file_format_str)) {
    @file_formats = split /,/, $file_format_str;
}

my @file_format_types;
if (defined($file_format_type_str)) {
    @file_format_types = split /,/, $file_format_type_str;
}

# Process the filter json string
my %json_filter;
if (defined($filter_json_str)) {
    my @tmp = split /,/, $filter_json_str;
    foreach my $term (@tmp) {
	my @tmp2 = split /=/, $term;
	if (exists($json_filter{$tmp2[0]})) {
	    push @{$json_filter{$tmp2[0]}}, $tmp2[1];
	} else {
	    $json_filter{$tmp2[0]} = [$tmp2[1]];
	}
    }
}

my %json_exclude;
if (defined($excl_json_str)) {
    my @tmp = split /,/, $excl_json_str;
    foreach my $term (@tmp) {
        my @tmp2 = split /=/, $term;
        if (exists($json_exclude{$tmp2[0]})) {
            push @{$json_exclude{$tmp2[0]}}, $tmp2[1];
        } else {
            $json_exclude{$tmp2[0]} = [$tmp2[1]];
        }
    }
}

# Set up the virtual browser
my $mech = WWW::Mechanize->new(
    ssl_opts => {
        verify_hostname => 0,
        # Quick and dirty hack to avoid errors about missing SSL certificates.
    },
    stack_depth => 0,
    autocheck => 0
    );

# Handle biosample-based queries with a preliminary query to get the
# biosample term name based on the search term.
my $biosample;
if ($by_biosample) {
    print STDERR "\nSearching ENCODE for a biosample matching $search_str...\n";
    $biosample = &get_biosample($mech, $search_str);
    if (!defined($biosample)) {
	die "Cannot find a biosample matching $search_str!\n";
    }
    $search_str = "*";
    $biosample = uri_escape($biosample);
} else {
    # Escape special characters in the search string
    $search_str = uri_escape($search_str);
}

# Get assembly params for certain types of queries from the URL params,
# $json_filter, and $json_exclude
my $assemblies = get_assembly_params($params, \%json_filter, \%json_exclude);

####################
# Stage 1: Build the query URL and run the primary query against the ENCODE
# Portal.


# Build the query URL from the command line args and parameters
my $URL = 'https://www.encodeproject.org/search/?searchTerm=';
$URL .= $search_str;
$URL .= "&type=";
$URL .= $rec_type;
$URL .= $params;

if (defined($biosample)) {
    $URL .= "&biosample_term_name=";
    $URL .= $biosample;
}

# &limit=all: Return all results, not just the first 25
# &frame=object: Return all non-null fields within the JSON data for results
# &format=json: Return search data in JSON format
$URL .= '&limit=all&frame=object&format=json';

print STDERR "Query URL: $URL\n";

# Retrieve the JSON data from the URL
my $json = &get_json($mech, $URL);

# Print status and result count to the terminal
print STDERR "${$json}{notification}: ${$json}{total} results found.\n\n";

if ($count_only) {
    print "${$json}{total}\n";
    exit 0;
}

# Catch requests that returned an error status from the Portal
if ( ${$json}{error_status} ) {
    exit 1;
}

#####################
# Stage Two: Find the files we need within the search results.
# Individual search results are within the embedded array @graph.

# Prepare the metadata file header
my @header;
if ($download || $file_list) {
    # File-centric header
    @header = ("file", "accession", "output_type", "biological_replicate",
               "technical_replicate", "assembly", "status", "dataset_type", "biosample_term_name",
               "biosample_type", "biosample_donor_accession", "assay_term_name", "target", "status",
               "date_released", "lab", "parent accession", "possible_controls",
               "documents", "href");
    if ($is_reads) {
	push @header, "run_type", "paired_end", "paired_with", "controlled_by";
    }
} elsif ($quality_data) {
    @header = ("experiment", "biosample_term_name",
               "biosample_type", "assay_term_name", "target", "file", "status",
               "assembly", "date_released", "lab", "files", "reads", "pctMapped",
	       "pctDup", "NSC", "RSC", "ccFragLen", "ccPhantomPeak",
	       "ccMin", "PBC1", "PBC2", "NRF");
} else {
    # Experiment-centric header
    @header = ("accession", "dataset_type", "biosample_term_name",
               "biosample_type", "assay_term_name", "target", "status",
               "date_released", "lab", "files", "possible_controls",
               "documents");
}

# Open a stream to the metadata file
my $DATA;
my $outfile = build_write_path($download_path, $out_root, ["metadata"]);
open $DATA, '>', $outfile;

# Print the metadata file header, if requested.
unless ($no_header) {
    &print_array(\@header, "\t", $DATA);
}


#my $mech_size = total_size($mech);
#print STDERR "Size of $mech: $mech_size\n";

my @downloads;  # To hold pointers to the files we will download
my @metadata;   # To hold metadata for the files we will download
my $n_files = 0;
my $n_results = 0;
my $ds = 1;

my %idx;
if ($n_rand > 0) {
    for (my $i = 0; $i < $n_rand; $i++) {
	my $m = int(rand(${$json}{total}));
	while (exists($idx{$m})) {
	    $m = int(rand(${$json}{total}));
	}
	$idx{$m} = 1;
    }
}


my $i = 0;
foreach my $row (@{${$json}{'@graph'}}) {

    if ($n_rand > 0) {
	if (!exists($idx{$i++})) {
	    next;
	}
    }
    
    if ($download || $file_list || $quality_data) {
	# If we are downloading data files, use a file-centric process and
	# metadata format. Quality data is the outlier here because it will
	# be compiled on an experiment-wise basis, but we still want to apply
	# the same file-wise filtering process.
	my $rec = $ds;
	if ($n_rand > 0) {
	    $rec = $i;
	}
	print STDERR "Processing files for result $rec...\n";
	$ds++;

	my @quality_recs; # Json objects from which to retrieve quality data

	my $has_rec = 0; # Does experiment have a useable record?
	my @files = @{$row->{files}};
	
	foreach my $file (@files) {
	    # Get json record for the file from ENCODE
	    my $file_json = get_file_json($file, $mech);
	    if (defined($file_json->{error_status}) &&
		$file_json->{error_status} == 1) {
		print STDERR "JSON retrieval error for $file: $file_json->{notification}.\n";
		next;
	    }

	    # If we are retrieving quality metrics, check to see if we have
	    # an alignment file that matches our assembly criteria.
	    if ($quality_data &&
		$file_json->{output_type} =~ m/alignments/ &&
		$file_json->{output_type} !~ m/unfiltered/) {
		push @quality_recs, {exp_json => $row,
				     file_json => $file_json};
	    }
	    
	    # Examine the JSON to see if the file matches our criteria
	    my $use_rec = 1; # Download all files by default

	    # If output_type does not match, reject the file
	    if (@output_types) {
		$use_rec = 0;
		for (my $i = 0; $i <= $#output_types ; $i++) {
		    if (screen_output_type($file_json, $output_types[$i])) {
			$use_rec = 1;
			last;
		    }
		}
	    }
	    
	    # If file_format does not match, reject the file
	    if (@file_formats && $use_rec) {
		for (my $i = 0; $i <= $#file_formats ; $i++) {
		    if (!exists(${$file_json}{file_format})) {
			$use_rec = 0;
			last;
		    } elsif (${$file_json}{file_format} eq $file_formats[$i]) {
			$use_rec = 1;
			last;
		    } else {
			if ($debug) {
			    print STDERR "file format ${$file_json}{file_format} does not match\n";
			}
			$use_rec = 0;
		    }
		}
	    }
		
	    # If file_format_type does not match, reject the file
	    if (@file_format_types && $use_rec) {
		for (my $i = 0; $i <= $#file_format_types ; $i++) {
		    if (!exists(${$file_json}{file_format_type})) {
			$use_rec = 0;
                        last;
		    } elsif (${$file_json}{file_format_type} eq $file_format_types[$i]) {
			$use_rec = 1;
			last;
		    } else {
			if ($debug) {
			    print STDERR "file format type ${$file_json}{file_format_type} does not match\n";
			}
			$use_rec = 0;
		    }
		}
	    }

	    # If we have other terms to check for in the json data, check them
	    if (%json_filter && $use_rec) {
		$use_rec = &filter_json(\%json_filter, $file_json, $match_all);
	    }

	    # If we have any terms to exclude, check for them
	    if (%json_exclude && $use_rec) {
		$use_rec = &exclude_json(\%json_exclude, $file_json, $match_all_excl);
            }

	    if ($use_rec && !$has_rec) {
		$has_rec = 1;
	    }

	    if ($use_rec && !$quality_data) {
		# Store a pointer to the file_json in our array
		if ($download) {
		    push @downloads, $file_json;
		}
		
		# Build a row for the metadata table
		my $controls_str = join ",", &nopaths($row->{possible_controls});
		my $documents_str = join ",", &nopaths($row->{documents});

		my $biorep = $file_json->{replicate}->{biological_replicate_number};
		if (!defined($biorep)) {
		    $biorep = $file_json->{biological_replicates}->[0];
		}
		my $trep = $file_json->{replicate}->{technical_replicate_number};
		if (!defined($trep)) {
		    $trep = $file_json->{technical_replicates}->[0];
		}

		my $url = 'http://www.encodeproject.org' . $file_json->{dataset} . '?format=json';
		#print STDERR "$url\n";
		my $pd_json = get_json($mech, $url);
		my $idx = 0;
		if ( defined($biorep) && defined($pd_json->{replicates}->[$biorep-1]) ) {
		    $idx = $biorep-1;
		}
		my $bs_acc = $pd_json->{replicates}->[$idx]->{library}->{biosample}->{donor}->{accession};

		my $biosample_ontology = get_biosample_ontology($mech, $row->{biosample_ontology});
		my @meta = (&get_filename($file_json),
			    $file_json->{accession},
			    $file_json->{output_type},
			    $biorep,
			    $trep,
			    $file_json->{assembly},
			    $file_json->{status},
			    $row->{dataset_type},
			    $biosample_ontology->{term_name},
			    $biosample_ontology->{classification},
			    $bs_acc,
			    $row->{assay_term_name},
			    &nopath($row->{target}),
			    $row->{status},
			    $row->{date_released}, 
			    &nopath($row->{lab}),
			    $row->{accession}, 
			    $controls_str,
			    $documents_str,
			    "https://www.encodeproject.org" . $file_json->{href},
		    );
		if ($is_reads) {
		    # Add in fastq read-specific values.
		    my $control_str = ".";
		    if (exists($file_json->{controlled_by})) {
			$control_str = join ",", nopaths($file_json->{controlled_by});
			if ($reads_and_controls) {
			    my $new_ctrls = parse_controls(\%controls,
							   $file_json->{controlled_by});
			    push @files, @{$new_ctrls};
			}
                    } else {
			if ($reads_and_controls && 
			    !exists($controls{"/files/" . $file_json->{accession}})) {
			    my $new_ctrls = find_controls($row->{possible_controls});
			    push @files, @{parse_controls(\%controls, $new_ctrls)};
			    # We need to prevent listing a file as its own control.
			    $control_str = join ",", nopaths($new_ctrls);
			}
		    }
		    
		    push @meta, $file_json->{run_type},
			$file_json->{paired_end},
			nopath($file_json->{paired_with}),
			$control_str;
		}
		
		# Correct any undefined metadata values
		for (my $i = 0; $i <= $#meta; $i++) {
		    if (!defined($meta[$i])) {
			$meta[$i] = '.';
		    }
		}

		# Print a metadata row
		&print_array(\@meta, "\t", $DATA);
		$n_results++
		
	    } else {
		# Should  not be necessary, but just to be sure there's no json
		# hanging around for files we're not using.
		undef($file_json);
	    }

	}
	
	if ($has_rec && $quality_data) {
	    my $res = get_quality_metrics(\@quality_recs, $params, \%json_filter, \%json_exclude, $mech, $assemblies);
	    #print STDERR "Line 699: success ", $res->{success}, ", meta ", $res->{meta},"\n";
	    if (!$res->{success}) {
		print STDERR "Metadata retrieval failed for $row->{accession}: $res->{message}\n";
		next;
	    }
	    foreach my $meta (@{$res->{res}}) {
		&print_array($meta, "\t", $DATA);
		$n_results++;
	    }
	}
	
    } else { # if (!$download)
	# Prepare a row for the tabular data...
	my $files_str = join ",", &nopaths($row->{files});
	my $controls_str = join ",", &nopaths($row->{possible_controls});
	my $documents_str = join ",", &nopaths($row->{documents});
	
	my @meta = ($row->{accession}, $row->{dataset_type},
		    $row->{biosample_term_name}, $row->{biosample_type},
		    $row->{assay_term_name}, &nopath($row->{target}),
		    $row->{status}, $row->{date_released},
		    &nopath($row->{lab}), $files_str, $controls_str,
		    $documents_str);
	for (my $i = 0; $i <= $#meta; $i++) {
	    if (!defined($meta[$i])) {
		$meta[$i] = '.';
	    }
	}
	print_array(\@meta, "\t", $DATA);
	$n_results++
    }
}

print STDERR "Wrote $n_results rows of metadata to \"$outfile\"....\n";

#$mech_size = total_size($mech);
#print STDERR "Size of $mech: $mech_size\n";

#######################
# Stage 3: Download data. Optionally save the file json.

# Download the files
if ($download && !$file_list) {
    print STDERR "Downloading $n_results files from ${$json}{total} records...\n";
    foreach my $file_json (@downloads) {
	&download_file($mech, $file_json, $download_path, $out_root, $check_exists, $use_wget);
	if ($save_json) {
	    &save_json($file_json, $download_path, $out_root);
	}
	#    $mech_size = total_size($mech);
	#    print STDERR "Size of $mech: $mech_size\n";
    }
}

close $DATA;

print STDERR "Done!\n";

######################
# Subroutines

sub get_json {
    # Retrieve JSON data from a URL and return it as a data structure.
    my ($mech, $url, $max_tries, $try) = @_;
    if (!defined($max_tries)) {
        $max_tries = 10;
    }
    if (!defined($try)) {
        $try = 1;
    }

    $mech->get($url);

    my $json;
    if ($mech->success()) {
        $json = decode_json($mech->content);
    } else {
        if ($try == $max_tries) {
            $json = decode_json('{ "notification": "Max tries reached. Aborting. This could be a problem with the ENCODE Portal, an empty result set, or a problem with your query. Please check your query for errors!", "total": 0, "error_status": 1 }');
        } else {
	    $json = get_json($mech, $url, $max_tries, ++$try);
        }
    }
    return $json;
}

sub save_json {
    # Save the JSON data for a file we have downloaded.
    my $json = $_[0];
    my $download_path = $_[1];
    my $out_root = $_[2];

    
    my $outfile = build_write_path($download_path, $out_root, [${$json}{accession}, 'json']);
    open my $OUT, '>', $outfile;
    print $OUT to_json($json, {pretty=>1});
    close $OUT;
}

sub build_write_path {
    # Build a filename and write path
    my ($download_path, $out_root, $args) = @_;
    # Last variable is an array reference for any additional file name parts
    # including the file extension. These will be '.' delimited.
    my $outfile;
    if (defined($out_root)) {
	$outfile = join '.', ($out_root, @{$args});
    } else {
	$outfile = join '.', @{$args};
    }
    if (defined($download_path)) {
	if ($download_path !~ m/\/$/) {
	    $download_path .= '/';
	}
	$outfile = $download_path . $outfile;
    }
    #print STDERR "$outfile";
    return $outfile;
}

sub download_file {
    # Given a file record, download the file from the ENCODE repository and
    # verify the MD5 checksum.
    my ($mech, $json, $download_path, $out_root, $check_exists, $use_wget) = @_;

    my @file_parts = split /\//, ${$json}{href};
    my $outfile = build_write_path($download_path, $out_root, [$file_parts[$#file_parts]]);

    if ($check_exists && -f $outfile) {
        # See if file already exists and move on if found
	print STDERR "\tFile exists: $outfile. Moving on...\n";
	return 1;
    }

    my $url;
    # Use S3 if possible.
    if (exists(${$json}{s3_uri})) {
	# Build an https download URL from the S3 URL
	my @parts = split /\//, ${$json}{s3_uri};
	$url = "https://" . $parts[2] . ".s3.amazonaws.com/" . join "/", @parts[3..$#parts];
    } else {
	$url = "https://www.encodeproject.org" . ${$json}{href};
    }
    print STDERR "Found a matching record at $url. Retrieving data...\n";

    if ($use_wget) {
	print STDERR "Retrieving file with wget...\n";
	# Retrieve with -q (quiet), -nc (no clobber), -c (continue download of partial files)
	`wget -q -nc -c $url`;
	if ($file_parts[$#file_parts] ne $outfile) {
	    `mv $file_parts[$#file_parts] $outfile`;
	}
	return 0;
    }

    $mech->get($url);
    
    # Make sure we got the file we were expecting by verifying the md5 checksum
    print STDERR "\tVerifying Checksum...\n";
    my $checksum = md5_hex($mech->content);
    if ($checksum ne ${$json}{md5sum}) {
	print STDERR "\tERROR: Checksums do not match. Skipping file!\n";
    } else {
	# Checksums match. Save file to the specified
	# path and file name.
	print STDERR "\tSaving file to $outfile...\n";
	$mech->save_content($outfile);
    }    
}

sub print_array {
    # A generic function to print an array to a file handle as a string with
    # the supplied delimiting.
    my ($array, $delim, $fh, $term, $nd_char) = @_;
    # input array, delimiter char, filehandle, terminator char, not-defined char                  

    if (!defined($delim)) {
        $delim = "\t";
    }
    if (!defined($term)) {
        $term = "\n";
    }

    if (defined($fh)) {
        select $fh;
    } else {
        select STDOUT;
    }

    my $i;
    for ($i = 0; $i < $#{$array}; $i++) {
        if (defined(${$array}[$i])) {
            print "${$array}[$i]$delim";
        } else {
            if (defined($nd_char)) {
                print "$nd_char$delim";
            } else {
                print STDERR "WARNING: Some fields in output array not defined with no default. Skipping!\n";
                next;
            }
        }
    }
    if (defined(${$array}[$#{$array}])) {
        print "${$array}[$#{$array}]$term";
    } else {
        if (defined($nd_char)) {
            print "$nd_char$term";
        } else {
            print STDERR "WARNING: Some fields in output array not defined with no default. Skipping!\n";
            print "$term";
        }
    }

    if (defined($fh)) {
        select STDOUT;
    }
    return 0;
}

sub nopaths {
    # A generic function to strip the path from file or directory names stored
    # in an array.
    my ($arr) = @_;

    my @out;
    for (my $i = 0; $i <= $#{$arr}; $i++) {
	my $s = &nopath($arr->[$i]);
	push @out, $s;
    }
    return @out;
}

sub nopath {
    # Strip the path from a file or directory name.
    my $str = $_[0];
    if (!defined($str)) {
	return ".";
    }
    $str =~ s/\/$//;
    $str =~ s/\/\S+\///;
    return $str;
}

sub get_biosample {
    my ($mech, $search_str) = @_;

    my $URL = 'http://www.encodeproject.org/search/?searchTerm=';
    $URL .= $search_str;
    $URL .= "&type=biosample";
    $URL .= '&limit=all&frame=object&format=json';
    
#    print STDERR "$URL\n";

    my $json = &get_json($mech, $URL);

    if ( ${$json}{error_status} ) {
	return undef;
    }

    my %terms_hash;
    foreach my $row (@{${$json}{'@graph'}}) {
	$terms_hash{${$row}{biosample_term_name}}++;
    }

    if (!%terms_hash) {
	return undef;
    }

    my @terms;
    foreach my $key (keys(%terms_hash)) {
	my @tmp = ($key, $terms_hash{$key});
	push @terms, \@tmp;
    }

    my $biosample;
    if ($#terms > 0) {
	print STDERR "More than one possible biosample found...\n\n";
	for (my $i = 0; $i <= $#terms; $i++) {
	    print STDERR $i+1, ": $terms[$i][0] ($terms[$i][1] occurences)\n";
	}
	print STDERR "\nEnter the number of the term you want to use [1] : ";
	my $idx = <STDIN>;
	if ($idx eq '') {
	    $idx = 1;
	}
	chomp $idx;
	$idx--;

	$biosample = $terms[$idx][0];

    } else {
	$biosample = $terms[0][0];
    }
    
    print STDERR "Using \"$biosample\" as biosample term name.\n\n";
    return $biosample;
}

sub get_biosample_ontology {
    my ($mech, $search_str) = @_;
    my $URL = 'http://www.encodeproject.org' . $search_str . '?format=json';

    my $json = &get_json($mech, $URL);

    if ( ${$json}{error_status} ) {
        return undef;
    }

    return $json;
}

sub get_filename {
    my ($file_json) = @_;
    if (exists($file_json->{s3_uri})) {
	my @tmp = split "/", $file_json->{s3_uri};
	return $tmp[$#tmp];
    }
    return &nopath($file_json->{href});
}

sub filter_json {
    my ($json_filter, $file_json, $match_all) = @_;

    foreach my $key (keys(%$json_filter)) {
	my @tmp = split /\./, $key;
	my $query_field = $tmp[$#tmp];
	my $query_terms = $json_filter->{$key};
	
	my $target_json = $file_json;
	if ($#tmp > 0) {
	    for (my $j = 0; $j < $#tmp; $j++) {
		if (exists($file_json->{$tmp[$j]})) {
		    $target_json = $file_json->{$tmp[$j]};
		}
	    }
	}
	
	if (!exists($target_json->{$query_field})) {
	    print STDERR "WARNING: Specified JSON attribute $key does not exist in file JSON. File will not be downloaded (accession = $file_json->{accession})!\n";
	    return 0;
	}

	for (my $i = 0; $i <= $#{$query_terms}; $i++) {
	    my $found = &filter_json_term($target_json, $query_field, $query_terms->[$i]);
	    if ($found && (!$match_all || $i == $#{$query_terms})) {
		return 1;
	    } elsif (!$found && $match_all) {
		return 0;
	    }
	}
    }
    return 0;
}

sub exclude_json {
    my ($json_exclude, $file_json, $match_all) = @_;

    foreach my $key (keys(%$json_exclude)) {
	my @tmp = split /\./, $key;
        my $query_field = $tmp[$#tmp];
        my $query_terms = $json_exclude->{$key};

        my $target_json = $file_json;
        if ($#tmp > 0) {
            for (my $j = 0; $j < $#tmp; $j++) {
                if (exists($file_json->{$tmp[$j]})) {
                    $target_json = $file_json->{$tmp[$j]};
                }
            }
        }

        if (!exists($target_json->{$query_field})) {
	    # If the given json key does not exist, accept the file.
            return 1;
        }

	my $n_found = 0;
        for (my $i = 0; $i <= $#{$query_terms}; $i++) {
            my $found = &filter_json_term($target_json, $query_field, $query_terms->[$i]);
            if ($found) {
		if (!$match_all) {
		    return 0;
		}
		# Handle the match_all case by comparing number
		# of terms matched to number of query terms.
		if (++$n_found == $#{$query_terms}) {
		    return 0;
		}
            }
        }
    }
    return 1;
}


sub filter_json_term {
    my ($target_json, $query_field, $query_term) = @_;

    if (ref $target_json->{$query_field} eq "ARRAY") {
	# If the target attribute is an array, look for any match                           
	# to the query term in the array.                                                   
	foreach my $attr (@{$target_json->{$query_field}}) {
	    if ($attr eq $query_term) {
		return 1;
	    }
	}
    } else {
	if ($target_json->{$query_field} eq $query_term) {
	    return 1;
	}
    }
    return 0;
}

sub get_file_json {
    # Retrieve json for a given file.
    my ($file, $mech) = @_;
    
    my $url = 'http://www.encodeproject.org' . $file . '?format=json';
    #print STDERR "$url\n";

    # Get the JSON from the url
    my $file_json = &get_json($mech, $url);
    return ($file_json);
}

sub get_quality_metrics {
    # Get quality metrics from an experiment
    my ($quality_recs, $params, $json_filter, $json_exclude, $mech, $assemblies) = @_;
    
    # Find the alignment file for the first/chosen assembly.
    my $success = 0;
    my $message = "No matching file(s)";
    my @res;

    foreach my $rec (@{$quality_recs}) {
	my $row = $rec->{exp_json};
	my $file_json = $rec->{file_json};
	#print STDERR $file_json->{accession}, "\n";

	if (!exists($file_json->{notes})) {
	    next;
	}
	
	# Check the assembly
	my $use_rec = 0;
	if ($#{$assemblies->{assembly}} > 0) {	    
	    for (my $i = 0; $i <= $#{$assemblies->{assembly}}; $i++) {
		if ($file_json->{assembly} eq $assemblies->{assembly}->[$i]) {
		    $use_rec = 1;
		    last;
		}
	    }
	    if (!$use_rec) {
		next;
	    }
	}
	if ($#{$assemblies->{not_assembly}} > 0) {
	    for	(my $i = 0; $i <= $#{$assemblies->{not_assembly}}; $i++) {
		if ($file_json->{assembly} eq $assemblies->{not_assembly}->[$i]) {
		    $use_rec = 0;
		    last;
		}
	    }
	    if (!$use_rec) {
		next;
            }
	}

	# Get the quality metrics. These are all available in the "notes" field.
	my $notes_json = decode_json($file_json->{notes});

	# Get biosample data (This is now a nested request in the experiment recs -- how annoying!)
	my $biosample_ontology = get_json($mech,
					  join("", "https://encodeproject.org", 
					       $row->{biosample_ontology}, 
					       "?format=json")
	    );
	if ($biosample_ontology->{error_status}) {
	    #print STDERR $biosample_ontology->{notification},"\n";
	    $biosample_ontology->{term_name} = $biosample_ontology->{classification} = ".";
	}

	my $pct_mapped = 0;
	if (exists($notes_json->{qc}->{qc})) {
	    $pct_mapped = $notes_json->{qc}->{qc}->{mapped}[0] / $notes_json->{qc}->{qc}->{in_total}[0]
	} else {
	    $pct_mapped = $notes_json->{qc}->{mapped}[0] / $notes_json->{qc}->{in_total}[0];
	}
	
        my @meta = ($row->{accession},
		    $biosample_ontology->{term_name},
		    $biosample_ontology->{classification},
		    $row->{assay_term_name}, &nopath($row->{target}),
		    $file_json->{accession},
		    $file_json->{status},
		    $file_json->{assembly},
		    $row->{date_released},
		    &nopath($row->{lab}),
		    $notes_json->{qc}->{qc}->{in_total}[0], 
		    $pct_mapped,
		    $notes_json->{qc}->{dup_qc}->{percent_duplication},
		    $notes_json->{qc}->{xcor_qc}->{phantomPeakCoef},
		    $notes_json->{qc}->{xcor_qc}->{relPhantomPeakCoef},
		    $notes_json->{qc}->{xcor_qc}->{corr_estFragLen},
		    $notes_json->{qc}->{xcor_qc}->{corr_phantomPeak},
		    $notes_json->{qc}->{xcor_qc}->{min_corr},
		    $notes_json->{qc}->{pbc_qc}->{PBC1},
		    $notes_json->{qc}->{pbc_qc}->{PBC2},
		    $notes_json->{qc}->{pbc_qc}->{NRF}
	    );
        for (my $i = 0; $i <= $#meta; $i++) {
            if (!defined($meta[$i])) {
                $meta[$i] = '.';
            }
        }
	push @res, \@meta;
	$success = 1;
    }        
    return {success => $success, message => $message, res => \@res};
}

sub get_assembly_params {
    # Get assembly specifications from URL params, %json_filter and %json_exclude.
    # Returns a hash of arrays.
    my ($params, $json_filter, $json_exclude) = @_;

    my @assembly;
    my @not_assembly;
    if ($params =~ m/assembly=(\w+)\&*/) {
        $assembly[0] = $1;
    } elsif (exists($json_filter->{assembly})) {
        for (my $i = 0; $i <= $#{$json_filter->{assembly}}; $i++) {
            push @assembly, $json_filter->{assembly}[$i];
	}
    }
    if (exists($json_exclude->{assembly})) {
        for (my $i = 0; $i <= $#{$json_exclude->{assembly}}; $i++) {
            push @not_assembly, $json_exclude->{assembly}[$i];
	}
    }
    
    return { assembly => \@assembly,
	     not_assembly => \@not_assembly };
}

sub parse_controls {
    my ($controls_found, $controls_this) = @_;
    my @ret;
    foreach my $control (@{$controls_this}) {
	# ENCODE is inconsistent in whether a trailing slash is
	# used in strings for controls and read files. Therefore,
	# we need the following hack to ensure nothing has a
	# trailing slash.
	$control =~ s/\/$//;
	if (!exists($controls_found->{$control})) {
	    $controls_found->{$control} = 1;
	    push @ret, $control;
	}
    }
    return \@ret;
}

sub find_controls {
    # Find control read files when they are not explicitly given in read file json
    my ($possible_controls) = @_;

    my @read_files;
    foreach my $pos_ctrl (@{$possible_controls}) {
	# It seems Possible Controls does not always include the /experiment/ prefix.
	if ($pos_ctrl !~ m/^\/experiments\//) {
	    $pos_ctrl = '/experiments/' . $pos_ctrl;
	}
	my $json = get_json($mech, 'http://www.encodeproject.org' . $pos_ctrl . '?format=json');
	foreach my $file_json (@{$json->{files}}) {
	    if (screen_output_type($file_json, "reads")) {
		push @read_files, '/files/' . $file_json->{accession};
	    }
	}
    }
    return \@read_files;
}

sub screen_output_type {
    # See if the output_type of the given record is what we want
    my ($file_json, $output_type) = @_;
    if (!exists($file_json->{output_type})) {
	return 0;
    } elsif ($file_json->{output_type} eq $output_type) {
	return 1;
    } else {
	if ($debug) {
	    print STDERR "Output type ${$file_json}{output_type} does not match specified type (", $output_type, ")\n";
	}
	return 0;
    }
}
