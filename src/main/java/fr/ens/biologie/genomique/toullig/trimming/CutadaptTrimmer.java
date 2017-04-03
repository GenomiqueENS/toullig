package fr.ens.biologie.genomique.toullig.trimming;

import com.google.common.base.Strings;
import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.Sequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import static fr.ens.biologie.genomique.eoulsan.bio.Alphabets.AMBIGUOUS_DNA_ALPHABET;
import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;
import static fr.ens.biologie.genomique.toullig.Utils.complement;

/**
 * Class to execute the Cutadapt trimmer. Created by birer on 27/03/17.
 */
public class CutadaptTrimmer implements Trimmer {

  private final Map<String, InformationRead> workTrimmingMap;
  private final File nameOutputFastq;
  private final File outputTrimLeftFastaFile;
  private final File outputTrimRightFastaFile;

  private final File infoTrimLeftFile;
  private final File infoTrimRightFile;

  private final String adaptorRetroTranscritpion;
  private final String adaptorStrandSwitching;
  private final Alphabet alphabet = AMBIGUOUS_DNA_ALPHABET;
  private final double errorRateCutadapt;

  private File fastaLeftOutlierFile;
  private File fastaRightOutlierFile;

  private FastaWriter leftFastaWriter;
  private FastaWriter rightFastaWriter;

  public CutadaptTrimmer(Map<String, InformationRead> workTrimmingMap,
      File nameOutputFastq, File outputTrimLeftFastaFile,
      File outputTrimRightFastaFile, String adaptorRetroTranscritpion,
      String adaptorStrandSwitching, double errorRateCutadapt,
      File fastaLeftOutlierFile, File fastaRightOutlierFile,
      File infoTrimLeftFile, File infoTrimRightFile) throws IOException {

    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = nameOutputFastq;
    this.outputTrimLeftFastaFile = outputTrimLeftFastaFile;
    this.outputTrimRightFastaFile = outputTrimRightFastaFile;
    this.adaptorRetroTranscritpion = adaptorRetroTranscritpion;
    this.adaptorStrandSwitching = adaptorStrandSwitching;
    this.errorRateCutadapt = errorRateCutadapt;
    this.infoTrimLeftFile = infoTrimLeftFile;
    this.infoTrimRightFile = infoTrimRightFile;
    this.fastaLeftOutlierFile = fastaLeftOutlierFile;
    this.fastaRightOutlierFile = fastaRightOutlierFile;

    this.leftFastaWriter = new FastaWriter(fastaLeftOutlierFile);
    this.rightFastaWriter = new FastaWriter(fastaRightOutlierFile);

  }

  //
  // Merge outlier
  //

  /**
   * Method of the class CutadaptTrimmer to merge the results of cutadapt with
   * the trim sequence.
   */
  public void mergeTrimOutlier() {

    HashMap<String, String> fastaLeftHash = new HashMap<>();
    HashMap<String, String> fastaRightHash = new HashMap<>();

    String shortestFastqSequence = "";
    int i = 0;
    int countWritten = 0;
    int countNull = 0;

    // Read of the left outlier fasta output by cutadapt

    try (FastaReader reader = new FastaReader(this.outputTrimLeftFastaFile)) {

      for (final Sequence read : reader) {
        fastaLeftHash.put(read.getName(), read.getSequence());
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    // Read of the right outlier fasta output by cutadapt

    try (FastaReader reader = new FastaReader(this.outputTrimRightFastaFile)) {

      for (final Sequence read : reader) {
        fastaRightHash.put(read.getName(), read.getSequence());
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    // Start merging
    System.out.println("Start of merging sequence !");

    try (BufferedWriter fastqTrimFile =
        new BufferedWriter(new FileWriter(this.nameOutputFastq))) {
      // merge of the outlier trim on the main sequence
      for (String id : this.workTrimmingMap.keySet()) {

        i++;
        String leftSequence;
        String rightSequence;
        InformationRead informationRead = this.workTrimmingMap.get(id);
        String cigar = informationRead.cigar;
        String quality = informationRead.quality;
        String sequence = informationRead.sequence;
        int leftLengthOutlier = informationRead.leftLengthOutlier;
        int rightLengthOutlier = informationRead.rightLengthOutlier;

        if (!cigar.equals("*")) {

          if (leftLengthOutlier < 0) {
            leftLengthOutlier = 0;
          }
          if (rightLengthOutlier < 0) {
            rightLengthOutlier = 0;
          }

          String mainSequenceWithoutOutlier = "";

          if ((leftLengthOutlier + rightLengthOutlier) < sequence.length()) {
            mainSequenceWithoutOutlier = sequence.substring(leftLengthOutlier,
                sequence.length() - rightLengthOutlier);
          }

          if (fastaLeftHash.containsKey(id)) {
            leftSequence = "";
          } else {
            leftSequence = fastaLeftHash.get(id);
          }
          if (fastaLeftHash.containsKey(id)) {
            rightSequence = "";
          } else {
            rightSequence = fastaRightHash.get(id);
          }

          // test if the sequence is null and put it empty if it's null
          rightSequence = Strings.nullToEmpty(rightSequence);

          leftSequence = Strings.nullToEmpty(leftSequence);

          int leftLength = leftLengthOutlier - leftSequence.length();
          int rightLength =
              (sequence.length() - rightLengthOutlier) + rightSequence.length();

          if (leftLengthOutlier - leftSequence.length() < 0) {
            leftLength = 0;
          }

          if (rightLength > quality.length()) {
            rightLength = quality.length();
          }

          // case with the addition of index (>0)
          if (leftLength > rightLength) {
            leftLength = rightLength;
          }

          String sequenceTrim = sequence.substring(leftLength, rightLength);

          String scoreTrim = quality.substring(leftLength, rightLength);

          if (scoreTrim.length() != sequenceTrim.length()) {
            System.out.println("problem :  "
                + scoreTrim.length() + "     " + sequenceTrim.length() + "   "
                + mainSequenceWithoutOutlier.length());
            System.out
                .println(leftLength + "     " + rightLength + "     " + id);
          }

          if (!sequenceTrim.isEmpty()) {

            ReadSequence fastq = new ReadSequence();
            fastq.setName(id);
            fastq.setSequence(sequenceTrim);
            fastq.setQuality(scoreTrim);
            FastqWriter fastqWriter = new FastqWriter(fastqTrimFile);
            fastqWriter.write(fastq);

            countWritten++;

            if (i == 1) {
              shortestFastqSequence = sequenceTrim;
            }
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

  //
  // Cutadapt
  //

  /**
   * Method of the class CutadaptTrimmer to execute cutadapt.
   * @param fastaOutlierFile, a fasta File with outlier to trim
   * @param strand, the cutadapt strand of the adaptor
   * @param infoTrimFile, path to the output log of cutadapt
   * @param pathOutputTrimFasta, path to the output of cutadapt
   */
  public void cutadaptTrimming(File fastaOutlierFile, String strand,
      File infoTrimFile, File pathOutputTrimFasta) {

    try {
      String reverseComplementAdaptorRT = strand
          + " reverse_complement_RT_adaptor="
          + reverseComplement(this.adaptorRetroTranscritpion, this.alphabet);
      String reverseComplementAdaptorSwithStrand = strand
          + " reverse_complement_Switch_Strand_RT_adaptor="
          + reverseComplement(this.adaptorStrandSwitching, this.alphabet);

      String complementAdaptorRT = strand
          + " complement_RT_adaptor="
          + complement(this.adaptorRetroTranscritpion, this.alphabet);
      String complementAdaptorSwithStrand = strand
          + " complement_Switch_Strand_RT_adaptor="
          + complement(this.adaptorStrandSwitching, this.alphabet);

      StringBuilder reverseAdaptorRTStringBuffer =
          new StringBuilder(this.adaptorRetroTranscritpion);
      StringBuilder reverseAdaptorSwithStrandStringBuffer =
          new StringBuilder(this.adaptorStrandSwitching);
      String reverseAdaptorRT = strand
          + " reverse_RT_adaptor="
          + reverseAdaptorRTStringBuffer.reverse().toString();
      String reverseAdaptorSwithStrand = strand
          + " reverse_Switch_Strand_RT_adaptor="
          + reverseAdaptorSwithStrandStringBuffer.reverse().toString();

      String adaptorRT =
          strand + " RT_adaptor=" + this.adaptorRetroTranscritpion;
      String adaptorStrandSwitching =
          strand + " Switch_Strand_RT_adaptor=" + this.adaptorStrandSwitching;

      String quiet = "";
      // String quiet="--quiet";

      ProcessBuilder pb = new ProcessBuilder("/bin/bash", "-c", "cutadapt "
          + adaptorRT + " " + adaptorStrandSwitching + " " + reverseAdaptorRT
          + " " + reverseAdaptorSwithStrand + " " + complementAdaptorRT + " "
          + complementAdaptorSwithStrand + " " + reverseComplementAdaptorRT
          + " " + reverseComplementAdaptorSwithStrand + " " + quiet
          + " --error-rate=" + this.errorRateCutadapt + " --info-file="
          + infoTrimFile.toString()
          + " --overlap=7 --times=8 --match-read-wildcards --format=fasta "
          + fastaOutlierFile.toString() + " > " + pathOutputTrimFasta);

      // System.out.println(pb.command());
      pb.redirectErrorStream(true);
      Process proc = pb.start(); // Start the process.

      getLogCutadapt(proc);

      proc.waitFor();

    } catch (IOException | InterruptedException e) {
      e.printStackTrace(); // or log it, or otherwise handle it
    }
  }

  /**
   * Method of the class CutadaptTrimmer to display log of cutadapt (delete the
   * option --quiet).
   * @param proc, a Processus
   * @throws IOException if an IO error occur
   */
  private void getLogCutadapt(Process proc) throws IOException {

    BufferedReader stdInput =
        new BufferedReader(new InputStreamReader(proc.getInputStream()));
    BufferedReader stdError =
        new BufferedReader(new InputStreamReader(proc.getErrorStream()));

    // read the output from the command
    System.out.println("Here is the standard output of the command:\n");
    String s;
    while ((s = stdInput.readLine()) != null) {
      System.out.println(s);
    }
    stdInput.close();

    // read any errors from the attempted command
    System.out.println("Here is the standard error of the command (if any):\n");
    while ((s = stdError.readLine()) != null) {
      System.out.println(s);
    }
    stdError.close();
  }

  /**
   * Method of the class CutadaptTrimmer to make some stats of the
   * cutadaptTrimming
   * @param infoTrimFile, a path to store stats in a file
   * @throws IOException if an IO error occur
   */
  public void statsLogCutadapt(File infoTrimFile) throws IOException {

    LocalReporter localReporterNumberTimesAdaptor = new LocalReporter();
    BufferedReader infoTrimBufferedReader =
        new BufferedReader(new FileReader(infoTrimFile));
    String line;
    String oldId = "";
    String stackConstructAdaptor = "";

    while ((line = infoTrimBufferedReader.readLine()) != null) {

      String[] part = line.split("\t");
      String id = part[0];
      String error = part[1];

      if (error.equals("-1")) {

        if (localReporterNumberTimesAdaptor.getCounterNames("Construction")
            .contains("no_adaptor_found")) {

          localReporterNumberTimesAdaptor.incrCounter("Construction",
              "no_adaptor_found", 1);
        } else {
          localReporterNumberTimesAdaptor.setCounter("Construction",
              "no_adaptor_found", 1);
        }
        continue;
      }

      String nameAdaptor = part[7];

      if (!id.equals(oldId)) {

        if (stackConstructAdaptor.equals("")) {
          if (localReporterNumberTimesAdaptor.getCounterNames("Construction")
              .contains(nameAdaptor)) {

            localReporterNumberTimesAdaptor.incrCounter("Construction",
                nameAdaptor, 1);
          } else {
            localReporterNumberTimesAdaptor.setCounter("Construction",
                nameAdaptor, 1);
          }
        } else {

          if (stackConstructAdaptor.contains(" + ")) {
            localReporterNumberTimesAdaptor.incrCounter("Construction",
                stackConstructAdaptor, 1);
          } else {
            localReporterNumberTimesAdaptor.incrCounter("Construction",
                stackConstructAdaptor, 1);
          }
          stackConstructAdaptor = "";
        }

      } else {
        if (stackConstructAdaptor.equals("")) {
          stackConstructAdaptor += nameAdaptor;
        } else {
          stackConstructAdaptor += " + " + nameAdaptor;
        }
      }
      oldId = id;
    }

    // Analyze counter group Construction
    for (String constructionAdaptor : localReporterNumberTimesAdaptor
        .getCounterNames("Construction")) {
      System.out.println(constructionAdaptor
          + " :  " + localReporterNumberTimesAdaptor
              .getCounterValue("Construction", constructionAdaptor));
    }
  }

  //
  // Main execution
  //

  /**
   * Method of the class TrimmomaticTrimmer to trimming with the Trimmer
   * interface.
   * @param leftLengthOutlier , the length of the left outlier
   * @param rightLengthOutlier , the length of the left outlier
   * @param sequence, the sequence of the read
   * @param id , the id of the read
   * @param quality , the quality of the read
   */
  public void preProcessTrimming(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality) {

    UtilsTrimming.writeOutliers(leftLengthOutlier, rightLengthOutlier, sequence,
        id, this.leftFastaWriter, this.rightFastaWriter);

  }

  /**
   * Method of the class TrimmomaticTrimmer to trimming with the Trimmer
   * interface
   */
  public void trimming() {

    System.out.println("Begin use cutadapt !");

    String strandLeft = "-g";
    String strandRight = "-a";

    // Cutadapt execution
    cutadaptTrimming(this.fastaLeftOutlierFile, strandLeft,
        this.infoTrimLeftFile, this.outputTrimLeftFastaFile);
    cutadaptTrimming(this.fastaRightOutlierFile, strandRight,
        this.infoTrimRightFile, this.outputTrimRightFastaFile);

    // Merge the output form cutadapt
    mergeTrimOutlier();

  }
}
