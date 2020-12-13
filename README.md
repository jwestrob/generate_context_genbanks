## generate_context_genbanks
Generating genbanks with annotations for use with Clinker.

---

# What does this script do?

This script takes protein fastas containing gene clusters of interest and their corresponding nucleotide fasta files, combines them with multiple possible annotation file types, and constructs genbank files ready to use with Clinker (https://github.com/gamcil/clinker).

# Why did you make a script to do that?

Because genbank files are a pain in the tuchus and figuring out the names of the subfields you can add annotations to is something you'd probably rather I do for you.

# What kind of annotation files can I use with this script?

You can use multiple annotation file types simultaneously, since there are multiple fields clinker parses ORF information from in a genbank file.

- pfamscan (http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/)

- kofamscan (https://github.com/takaram/kofam_scan)

- GOOSOS (https://github.com/jwestrob/GOOSOS)

- Custom annotation files (see below example and the examples folder in this repository)

# You talked about examples but I don't see any! I want to speak with your manager right now young man! Go get them. I'll wait.

Chill Karen. I need to make a gene cluster fasta and generate associated annotations for a public-domain genome that's not part of one of my ongoing projects, so give me a hot second OK?

## I keep trying to use kofamscan but your script gets mad at me or messes up in some miscellaneous way.

I'm not surprised; in the Banfield lab we've developed a processing pipeline that takes the output of KOFAMscan and selects only the above-threshold hits. If you're trying to use this script, want to include KOFAM info, and aren't in the Banfield lab, hit me up.

## What the heck? Custom annotation files?

Yes indeed. Please make this a tab separated, two column file with the first column being the protein IDs and the second the custom annotation you wish to display. See this example:

```
example_genome_id|example_orf_id_1  example_annotation_1
example_genome_id|example_orf_id_2  example_annotation_2
```

You might notice that I have the genome ID in the protein ID, separated with a `|` character. As of writing this (the first day of writing this script), I haven't explored a good number of edge cases, so if you don't have those genome IDs appended, there may be some key error issues. I'm going to make sure the script is robust to this in the future, but for now, caution when making custom annotations if you do not have the genome ID included.

---

# Wow, you're so cool Jacob! But your script sucks and won't run.

Yeah I haven't been working on this one for very long, and I'm sure there are a bunch of ways to break it since I'm dealing with numerous different file formats. If you experience problems send me a slack (if you're in the Banfield lab) or raise an issue here and I'll fix it PDQ.

