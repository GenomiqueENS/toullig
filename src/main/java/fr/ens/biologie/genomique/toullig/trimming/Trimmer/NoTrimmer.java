package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

/**
 * Created by birer on 14/04/17.
 */
public class NoTrimmer implements Trimmer {

    private final Map<String, InformationRead> workTrimmingMap;
    private final File nameOutputFastq;


    public NoTrimmer(Map<String, InformationRead> workTrimmingMap,
                     File outputFastqFile){

        this.workTrimmingMap = workTrimmingMap;
        this.nameOutputFastq = outputFastqFile;
    }


    private void writeSequenceWithNoOutlier(){

        int i=0;
        int countWritten=0;
        int countNull=0;
        String shortestFastqSequence="";

        // test if the output left fasta file is correctly open
        try (BufferedWriter fastqTrimFile =
                     new BufferedWriter(new FileWriter(this.nameOutputFastq))) {

            // get id for each red on the work map
            for (String id : this.workTrimmingMap.keySet()) {

                i++;

                // get information for the read
                InformationRead informationRead = this.workTrimmingMap.get(id);
                String cigar = informationRead.cigar;
                String quality = informationRead.quality;
                String sequence = informationRead.sequence;
                int leftLengthOutlier = informationRead.leftLengthOutlier;
                int rightLengthOutlier = informationRead.rightLengthOutlier;

                // test if the read is unmapped
                if (!cigar.equals("*")) {

                    // test if the length of the left outlier is negatif
                    if (leftLengthOutlier < 0) {
                        leftLengthOutlier = 0;
                    }

                    // test if the length of the left outlier is negatif
                    if (rightLengthOutlier < 0) {
                        rightLengthOutlier = 0;
                    }

                    String mainSequenceWithoutOutlier = "";

                    String sequenceTrim ="";
                    String qualityTrim = "";

                    // test if the rightlengthsequence and the leftlengthsequence are overlap
                    if((sequence.length()-rightLengthOutlier)>leftLengthOutlier){
                        sequenceTrim = sequence.substring(leftLengthOutlier, sequence.length()-rightLengthOutlier);
                        qualityTrim = quality.substring(leftLengthOutlier, quality.length()-rightLengthOutlier);
                    }


                    // test if the quality length is differerent to the sequence length
                    if (qualityTrim.length() != sequenceTrim.length()) {
                        System.out.println("problem :  "
                                + qualityTrim.length() + "     " + sequenceTrim.length() + "   "
                                + mainSequenceWithoutOutlier.length());
                        System.out
                                .println(leftLengthOutlier + "     " + rightLengthOutlier + "     " + id);
                    }

                    // test if the sequence trimmed is empty
                    if (!sequenceTrim.isEmpty()) {

                        ReadSequence fastq = new ReadSequence();
                        fastq.setName(id);
                        fastq.setSequence(sequenceTrim);
                        fastq.setQuality(qualityTrim);

                        // write the trimmed read
                        FastqWriter fastqWriter = new FastqWriter(fastqTrimFile);
                        fastqWriter.write(fastq);

                        countWritten++;

                        // test for the first loop
                        if (i == 1) {
                            shortestFastqSequence = sequenceTrim;
                        }

                        // test if a new shortest read is comput in this loop
                        if (shortestFastqSequence.length() >= sequenceTrim.length()) {
                            shortestFastqSequence = sequenceTrim;
                        }

                    } else {
                        countNull++;
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(
                "The shortest read size is: " + shortestFastqSequence.length());
        System.out.println("Number of trim read write: " + countWritten);
        System.out.println("Number of trim read null: " + countNull);

    }


    // No preprocessTrimming are require
    @Override
    public void preProcessTrimming(int leftLengthOutlier, int rightLengthOutlier, String sequence, String id, String quality) {

    }

    /**
     * Method of the class NoTrimmer to trimming with the Trimmer
     * interface
     */
    public void trimming() {

        writeSequenceWithNoOutlier();

    }

}