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
my $output_type_str;
my $file_format_str;
my $file_format_type_str;
my $download_path;
my $out_root;
my $save_json = 0;
my $rec_type = "experiment";
my $no_header = 0;
my $count_only = 0;
my $by_biosample = 0;
my $check_exists = 0;
my $use_wget = 0;
my $filter_json_str;
my $debug = 0;
my $n_rand = 0;

GetOptions (
    "help" => \$help,
    "download" => \$download,
    "file-list" => \$file_list,
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
    "random=i" => \$n_rand,
    "debug" => \$debug
    );

# Check for proper usage and help option and exit with usage message as needed.
if ($help || $#ARGV+1 < 2) {
    die "$usage\n";
}

# Check for other options and do some sanity checks...
if (defined($output_type_str) && !($download || $file_list) ) {
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
if ($rec_type ne "experiment") {
    print STDERR "\nWarning: --rec_type other than \"experiment\" has not been well tested and may produces unpredictable results!\n";
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
if (defined($filter_json_str) && !($download || $file_list)) {
    print STDERR "\nWarning: --filter_json_str has no effect without --download or --file-list.\n";
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
	$json_filter{$tmp2[0]} = $tmp2[1];
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

####################
# Stage 1: Build the query URL and run the primary query against the ENCODE
# Portal.

my $URL;

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

# Build the query URL from the command line args and parameters
$URL = 'https://www.encodeproject.org/search/?searchTerm=';
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

    # Each row of the array contains a hash reference to a search result.
    # Make a copy of the hash for convenience.
    my %result = %{$row};

    if ($download || $file_list) {
	# If we are downloading data files, use a file-centric process and
	# metadata format.
	my $rec = $ds;
	if ($n_rand > 0) {
	    $rec = $i;
	}
	print STDERR "Processing files for result $rec...\n";
	$ds++;
	
	my @files = @{$result{files}};
	foreach my $file (@files) {
	    
	    # Because metadata for individual files is stored in separate
	    # records, retrieve these as separate queries against the
	    # database.
	    my $url = 'http://www.encodeproject.org/' . $file . '?format=json';
	    #print STDERR "$url\n";
	    
	    # Get the JSON from the url
	    my $file_json = &get_json($mech, $url);
	    if (defined($file_json->{error_status}) &&
		$file_json->{error_status} == 1) {
		print STDERR "JSON retrieval error for $url: $file_json->{notification}.\n";
		next;
	    }
	    
	    # Examine the JSON to see if the file matches our criteria
	    my $use_rec = 1; # Download all files by default
		
	    # If output_type does not match, reject the file
	    if (@output_types) {
		for (my $i = 0; $i <= $#output_types ; $i++) {
		    if (!exists(${$file_json}{output_type})) {
			$use_rec = 0;
			last;
		    } elsif (${$file_json}{output_type} eq $output_types[$i]) {
			$use_rec = 1;
			last;
		    } else {
			if ($debug) {
			    print STDERR "output type ${$file_json}{output_type} does not match\n";
			}
			$use_rec = 0;
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
		foreach my $key (keys(%json_filter)) {
		    my @tmp = split /\./, $key;
		    
		    my $tmp_json = $file_json;
		    if ($#tmp > 0) {
			for (my $j = 0; $j < $#tmp; $j++) {
			    if (exists($file_json->{$tmp[$j]})) {
				$tmp_json = $file_json->{$tmp[$j]};
			    }
			}
		    }
		    
		    if (ref $tmp_json->{$tmp[$#tmp]} eq "ARRAY") {
			# If the target attribute is an array, look for any match
			# to the query term in the array.
			$use_rec = 0;
			foreach my $attr (@{$tmp_json->{$tmp[$#tmp]}}) {
			    if ($debug) {
				print STDERR "$attr, $json_filter{$key}\n";
			    }
			    if ($attr eq $json_filter{$key}) {
				$use_rec = 1;
			    } 
			}
		    } else {			
			if ($tmp_json->{$tmp[$#tmp]} ne $json_filter{$key}) {
			    $use_rec = 0;
			}
		    }	   
		    if (!exists($tmp_json->{$tmp[$#tmp]})) {
			$use_rec = 0;
			print STDERR "WARNING: Specified JSON attribute $key does not exist in file JSON. File will not be downloaded (accession = $file)!\n";
		    }
		}
	    }
	    
	    if ($use_rec) {
		# Store a pointer to the file_json in our array
		if ($download) {
		    push @downloads, $file_json;
		}
		
		# Build a row for the metadata table
		my $controls_str = join ",", &nopaths($result{possible_controls});
		my $documents_str = join ",", &nopaths($result{documents});

		my $biorep = $file_json->{replicate}->{biological_replicate_number};
		if (!defined($biorep)) {
		    $biorep = $file_json->{biological_replicates}->[0];
		}
		my $trep = $file_json->{replicate}->{technical_replicate_number};
		if (!defined($trep)) {
		    $trep = $file_json->{technical_replicates}->[0];
		}

		my $url = 'http://www.encodeproject.org' . $file_json->{dataset} . '?format=JSON';
		#print STDERR "$url\n";
		my $pd_json = get_json($mech, $url);
		my $idx = 0;
		if ( defined($biorep) && defined($pd_json->{replicates}->[$biorep-1]) ) {
		    $idx = $biorep-1;
		}
		my $bs_acc = $pd_json->{replicates}->[$idx]->{library}->{biosample}->{donor}->{accession};
		    
		my @row = (&nopath($file_json->{href}),
			   $file_json->{accession}, $file_json->{output_type},
			   $biorep,
			   $trep,
			   $file_json->{assembly},
			   $file_json->{status},
			   $result{dataset_type}, $result{biosample_term_name},
			   $result{biosample_type},
			   $bs_acc,
			   $result{assay_term_name},
			   &nopath($result{target}), $result{status},
			   $result{date_released}, &nopath($result{lab}),
			   $result{accession}, $controls_str, $documents_str,
			   "https://www.encodeproject.org" . $file_json->{href},
		    );

		# Correct any undefined metadata values
		for (my $i = 0; $i <= $#row; $i++) {
		    if (!defined($row[$i])) {
			$row[$i] = '.';
		    }
		}

		# Print a metadata row
		&print_array(\@row, "\t", $DATA);

		$n_results++
		
	    } else {
		# Should  not be necessary, but just to be sure there's no json
		# hanging around for files we're not using.
		undef($file_json);
	    }
	}
	
    } else { # if (!$download)
	# Prepare a row for the tabular data...
	my $files_str = join ",", &nopaths($result{files});
	my $controls_str = join ",", &nopaths($result{possible_controls});
	my $documents_str = join ",", &nopaths($result{documents});
	
	my @row = ($result{accession}, $result{dataset_type},
		   $result{biosample_term_name}, $result{biosample_type},
		   $result{assay_term_name}, &nopath($result{target}),
		   $result{status}, $result{date_released},
		   &nopath($result{lab}), $files_str, $controls_str,
		   $documents_str);
	for (my $i = 0; $i <= $#row; $i++) {
	    if (!defined($row[$i])) {
		$row[$i] = '.';
	    }
	}
	print_array(\@row, "\t", $DATA);
	$n_results++
    }
}

print STDERR "Wrote $n_results rows of metadata to \"$outfile\"....\n";

#$mech_size = total_size($mech);
#print STDERR "Size of $mech: $mech_size\n";

#######################
# Stage 3: Download data. Optionally save the file json.

# Download the files
if ($download) {
    print STDERR "Downloading $n_results files from ${$json}{total} records...\n";
}
foreach my $file_json (@downloads) {
    &download_file($mech, $file_json, $download_path, $out_root, $check_exists, $use_wget);
    if ($save_json) {
	&save_json($file_json, $download_path, $out_root);
    }
#    $mech_size = total_size($mech);
#    print STDERR "Size of $mech: $mech_size\n";
}

close $DATA;

print STDERR "Done!\n";

######################
# Subroutines

sub get_json {
    # Retrieve JSON data from a URL and return it as a data structure.
    my $mech = $_[0];
    my $url = $_[1];

    my $max_tries = 10; # TO-DO: make this configurable

    my $status = eval {
	$mech->get($url);
	1
    };

    my $content;
    if (!$status) {
	# For some reason, the portal sometimes returns an error when there
	# are no results. In these cases, use a dummy json to notify the user.
	return decode_json('{ "notification": "ENCODE Portal returned error status: no results or bad request. Please check your query for errors!", "total": 0, "error_status": 1 }');
    } else {	
	$content = $mech->content;
    }

    my $json;
    eval { $json = decode_json($content); };
    
    # Retry if there is an error.
    my $tries = 2;
    while ($@ && $tries <= $max_tries) {
	if ($debug) {
	    print STDERR "Retrying json retrieval: attempt $tries of $max_tries...\n";
	}
	$tries++;
	$content = $mech->content;
	eval { $json = decode_json($content); };
    }
    if ($@) {
	return decode_json('{ "notification": "Max tries reached. Aborting.", "total": 0, "error_status": 1 }');
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
    my $mech = $_[0];
    my $json = $_[1];
    my $download_path = $_[2];
    my $out_root = $_[3];
    my $check_exists = $_[4];
    my $use_wget = $_[5];


    my @file_parts = split /\//, ${$json}{href};
    my $outfile = build_write_path($download_path, $out_root, [$file_parts[$#file_parts]]);

    if ($check_exists && -f $outfile) {
        # See if file already exists and move on if found
	print STDERR "\tFile exists: $outfile. Moving on...\n";
	return 1;
    }

    my $url = "https://www.encodeproject.org" . ${$json}{href};
    print STDERR "Found a matching record at $url. Retrieving data...\n";

    if ($use_wget) {
	print STDERR "Retrieving file with wget...\n";
	`wget $url`;
	`mv $file_parts[$#file_parts] $outfile`;
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
    my $array = $_[0];
    my $delim = $_[1];
    my $fh = $_[2];

    my $out_str = join $delim, @{$array};
    print $fh "$out_str\n";
}

sub nopaths {
    # A generic function to strip the path from file or directory names stored
    # in an array.
    my @array = @{$_[0]};

    my @out;
    if ($#array < 0) {
	push @out, ".";
    } else {
	foreach my $str (@array) {
	    my $s = &nopath($str);
	    push @out, $s;
	}
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
