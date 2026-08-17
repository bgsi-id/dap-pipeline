#!/usr/bin/env perl
use strict;
use warnings;

# Add deterministic, typed rescue fields to a normalized site VCF.  The input is
# plain VCF on stdin and the output is plain VCF on stdout; bcftools owns BGZF and
# indexing in the Nextflow process.

sub info_map {
    my ($raw) = @_;
    my %result;
    return %result if !defined($raw) || $raw eq '.';
    for my $item (split /;/, $raw) {
        my ($key, $value) = split /=/, $item, 2;
        $result{$key} = defined($value) ? $value : '1';
    }
    return %result;
}

sub spliceai_max {
    my ($raw) = @_;
    return undef if !defined($raw) || $raw eq '.';
    my $maximum;
    for my $entry (split /,/, $raw) {
        my @parts = split /\|/, $entry, -1;
        next if @parts < 6;
        for my $index (2 .. 5) {
            next if $parts[$index] eq '' || $parts[$index] eq '.';
            next if $parts[$index] !~ /^(?:0(?:\.\d+)?|1(?:\.0+)?)$/;
            my $value = 0 + $parts[$index];
            $maximum = $value if !defined($maximum) || $value > $maximum;
        }
    }
    return $maximum;
}

sub approved_clinvar_rescue {
    my ($significance, $conflicts) = @_;
    for my $raw (grep { defined($_) && $_ ne '.' } ($significance, $conflicts)) {
        for my $token (split /[,\|\/]/, $raw) {
            $token =~ s/\([^)]*\)//g;
            $token =~ s/^\s+|\s+$//g;
            $token = lc($token);
            $token =~ s/[ -]+/_/g;
            return 1 if $token eq 'pathogenic' || $token eq 'likely_pathogenic';
        }
    }
    return 0;
}

while (my $line = <STDIN>) {
    if ($line =~ /^#CHROM\b/) {
        print qq{##INFO=<ID=DAP_SPLICEAI_DS_MAX,Number=1,Type=Float,Description="Maximum of SpliceAI DS_AG, DS_AL, DS_DG and DS_DL">\n};
        print qq{##INFO=<ID=DAP_CLINVAR_RESCUE,Number=0,Type=Flag,Description="Approved ClinVar P/LP detail-annotation rescue">\n};
        print $line;
        next;
    }
    if ($line =~ /^#/) {
        print $line;
        next;
    }

    chomp $line;
    my @fields = split /\t/, $line, -1;
    die "malformed VCF row\n" if @fields < 8;
    my %info = info_map($fields[7]);
    my @extra;
    my $maximum = spliceai_max($info{SpliceAI});
    push @extra, sprintf('DAP_SPLICEAI_DS_MAX=%.6g', $maximum) if defined($maximum);
    push @extra, 'DAP_CLINVAR_RESCUE'
        if approved_clinvar_rescue($info{CLNSIG}, $info{CLNSIGCONF});
    if (@extra) {
        $fields[7] = $fields[7] eq '.' ? join(';', @extra) : join(';', $fields[7], @extra);
    }
    print join("\t", @fields), "\n";
}
