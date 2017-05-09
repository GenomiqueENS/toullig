package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import com.google.common.base.Strings;
import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.Sequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import fr.ens.biologie.genomique.toullig.Utils;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;
import fr.ens.biologie.genomique.toullig.trimming.UtilsTrimming;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;
import static fr.ens.biologie.genomique.eoulsan.bio.Alphabets.AMBIGUOUS_DNA_ALPHABET;
import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;
import static fr.ens.biologie.genomique.toullig.Utils.complement;

/**
 * Class to execute the Cutadapt trimmer. Created by birer on 27/03/17.
 * @author Aurelien Birer
 */
class CutadaptTrimmer implements Trimmer {

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

  private final File fastaLeftOutlierFile;
  private final File fastaRightOutlierFile;

  /**
   * Constructor of the CutadaptTrimmer class
   * @param workTrimmingMap, a Map of working information
   * @param outputFastqFile, a File of the output fastq
   * @param outputTrimLeftFastaFile, a File of the left outlier fasta
   * @param outputTrimRightFastaFile, a File of the right outlier fasta
   * @param adaptorRetroTranscritpion, a string of the Retro-Transcription
   *          adaptor
   * @param adaptorStrandSwitching, a string of the Strand-Switching adaptor
   * @param errorRateCutadapt, a int of the error rate for Cutadapt
   * @param fastaLeftOutlierFile, a File of the left outlier output fasta
   * @param fastaRightOutlierFile, a File of the right outlier output fasta
   * @param infoTrimLeftFile, a File of the left information give by cutadapt
   * @param infoTrimRightFile, a File of the right information give by cutadapt
   * @throws IOException, an IOException if an error occur
   */
  CutadaptTrimmer(Map<String, InformationRead> workTrimmingMap,
      File outputFastqFile, File outputTrimLeftFastaFile,
      File outputTrimRightFastaFile, String adaptorRetroTranscritpion,
      String adaptorStrandSwitching, double errorRateCutadapt,
      File fastaLeftOutlierFile, File fastaRightOutlierFile,
      File infoTrimLeftFile, File infoTrimRightFile) throws IOException {

    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = outputFastqFile;
    this.outputTrimLeftFastaFile = outputTrimLeftFastaFile;
    this.outputTrimRightFastaFile = outputTrimRightFastaFile;
    this.adaptorRetroTranscritpion = adaptorRetroTranscritpion;
    this.adaptorStrandSwitching = adaptorStrandSwitching;
    this.errorRateCutadapt = errorRateCutadapt;
    this.infoTrimLeftFile = infoTrimLeftFile;
    this.infoTrimRightFile = infoTrimRightFile;
    this.fastaLeftOutlierFile = fastaLeftOutlierFile;
    this.fastaRightOutlierFile = fastaRightOutlierFile;

  }

  /**
   * Method of the CutadaptTrimmer class to read a Cutadapt fasta output.
   * @param outputTrimFastaFile, a fasta output of cutadapt
   * @param fastaHash, a HashMap
   * @return , a HashMap
   */
  private HashMap<String, String> FastaReaderForCutadapt(
      File outputTrimFastaFile, HashMap<String, String> fastaHash) {

    // test if the output left fasta file is correctly open
    try (FastaReader reader = new FastaReader(outputTrimFastaFile)) {

      // get each read in the output
      for (final Sequence read : reader) {
        fastaHash.put(read.getName(), read.getSequence());
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    return fastaHash;
  }

  /**
   * Method of the class CutadapatTrimmer to get the outlier sequence.
   * @param fastaLeftHash, a fastaLeftHash with all outliers sequence
   * @param id, the id to the interest sequence
   * @return , the sequence of the outlier trim
   */
  public String lengthSequence(HashMap<String, String> fastaLeftHash,
      String id) {

    String lengthSequence;

    // test if the local left hash contains this id read
    if (fastaLeftHash.containsKey(id)) {
      lengthSequence = "";
    } else {
      lengthSequence = fastaLeftHash.get(id);
    }

    // test if the sequence is null and put it empty if it's null for the
    // left outlier
    lengthSequence = Strings.nullToEmpty(lengthSequence);

    return lengthSequence;
  }

  //
  // Merge outlier
  //

  /**
   * Method of the class CutadaptTrimmer to merge the results of cutadapt with
   * the trim sequence.
   */
  private void mergeTrimOutlier() {

    HashMap<String, String> fastaLeftHash = new HashMap<>();
    HashMap<String, String> fastaRightHash = new HashMap<>();

    String shortestFastqSequence = "";
    int i = 0;
    int countWritten = 0;
    int countNull = 0;

    //
    // Read of the outliers fastas outputs by cutadapt
    //

    fastaLeftHash =
        FastaReaderForCutadapt(this.outputTrimLeftFastaFile, fastaLeftHash);
    fastaRightHash =
        FastaReaderForCutadapt(this.outputTrimRightFastaFile, fastaRightHash);

    //
    // Start merging
    //

    getLogger().info("Start of merging sequence !");

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

          // test if the trim sequence and strictly inferior to the length
          // of the sequence
          if ((leftLengthOutlier + rightLengthOutlier) < sequence.length()) {

            mainSequenceWithoutOutlier = sequence.substring(leftLengthOutlier,
                sequence.length() - rightLengthOutlier);

          }

          // get the left sequence outlier trimmed
          String leftSequence = lengthSequence(fastaLeftHash, id);

          // get the right sequence outlier trimmed
          String rightSequence = lengthSequence(fastaRightHash, id);

          //
          // get the length of the trim outlier sequence
          //

          int leftLength = leftLengthOutlier - leftSequence.length();
          int rightLength =
              (sequence.length() - rightLengthOutlier) + rightSequence.length();

          // test if the left trim outlier sequence is superior to the left
          // outlier sequence
          if (leftLengthOutlier - leftSequence.length() < 0) {
            leftLength = 0;
          }

          // test if the right length trimmed sequence is superior to the
          // quality length
          if (rightLength > quality.length()) {
            rightLength = quality.length();
          }

          // case with the addition of index (>0)
          if (leftLength > rightLength) {
            leftLength = rightLength;
          }

          //
          // Get the sequence and quality trimmed
          //

          String sequenceTrim = sequence.substring(leftLength, rightLength);

          String qualityTrim = quality.substring(leftLength, rightLength);

          // test if the quality length is differerent to the sequence length
          if (qualityTrim.length() != sequenceTrim.length()) {
            System.out.println("problem :  "
                + qualityTrim.length() + "     " + sequenceTrim.length() + "   "
                + mainSequenceWithoutOutlier.length());
            System.out
                .println(leftLength + "     " + rightLength + "     " + id);
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

    getLogger()
        .info("The shortest read size is: " + shortestFastqSequence.length());
    getLogger().info("Number of trim read write: " + countWritten);
    getLogger().info("Number of trim read null: " + countNull);
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
  private void cutadaptTrimming(File fastaOutlierFile, String strand,
      File infoTrimFile, File pathOutputTrimFasta) {

    // test if the follow code execute correctly
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

      String reverseAdaptorRT = strand
          + " reverse_RT_adaptor="
          + Utils.reverse(this.adaptorRetroTranscritpion);

      String reverseAdaptorSwithStrand = strand
          + " reverse_Switch_Strand_RT_adaptor="
          + Utils.reverse(this.adaptorStrandSwitching);

      String adaptorRT =
          strand + " RT_adaptor=" + this.adaptorRetroTranscritpion;

      String adaptorStrandSwitching =
          strand + " Switch_Strand_RT_adaptor=" + this.adaptorStrandSwitching;

      // if a problem occur for the cutadapt display
      // String quiet="--quiet";
      String quiet = "";

      // get the processus builder for launching cutadapt in the shell
      ProcessBuilder pb = new ProcessBuilder("/bin/bash", "-c", "cutadapt "
          + adaptorRT + " " + adaptorStrandSwitching + " " + reverseAdaptorRT
          + " " + reverseAdaptorSwithStrand + " " + complementAdaptorRT + " "
          + complementAdaptorSwithStrand + " " + reverseComplementAdaptorRT
          + " " + reverseComplementAdaptorSwithStrand + " " + quiet
          + " --error-rate=" + this.errorRateCutadapt + " --info-file="
          + infoTrimFile.toString()
          + " --overlap=8 --times=10 --match-read-wildcards --format=fasta "
          + fastaOutlierFile.toString() + " > " + pathOutputTrimFasta);

      // System.out.println(pb.command());

      // get error enable for the processus
      pb.redirectErrorStream(true);

      // process the processus
      Process proc = pb.start(); // Start the process.

      // get log of the processus
      getLogCutadapt(proc);

      // wait to the end of the processus
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
  void statsLogCutadapt(File infoTrimFile, String header) throws IOException {

    // create new localReporter object
    LocalReporter localReporterNumberTimesAdaptor = new LocalReporter();

    // open to read the info Trim file of cutadapt
    BufferedReader infoTrimBufferedReader =
        new BufferedReader(new FileReader(infoTrimFile));
    String line;
    String oldId = "";
    String stackConstructAdaptor = "";

    // read the info trim file of cutadapt
    while ((line = infoTrimBufferedReader.readLine()) != null) {

      String[] part = line.split("\t");
      String id = part[0];
      String error = part[1];

      // test if no adaptor is found for this read by cutadapt
      if (error.equals("-1")) {

        // test if this status already exist on the localReporter object
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

      // test if the id uis the same of the last loop (technique to have multi
      // adaptor found)
      if (!id.equals(oldId)) {

        // test if the stack construct adaptor is empty
        if (stackConstructAdaptor.equals("")) {

          // test if the localReporterNumberTimesAdaptor have already this
          // adaptor construction
          if (localReporterNumberTimesAdaptor.getCounterNames("Construction")
              .contains(nameAdaptor)) {

            localReporterNumberTimesAdaptor.incrCounter("Construction",
                nameAdaptor, 1);
          } else {
            localReporterNumberTimesAdaptor.setCounter("Construction",
                nameAdaptor, 1);
          }
        } else {

          // test if the localReporterNumberTimesAdaptor have a complex adaptor
          // construction (more than 2 adaptor)
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

        // test if the stack construction adaptor is empty
        if (stackConstructAdaptor.equals("")) {
          stackConstructAdaptor += nameAdaptor;
        } else {
          stackConstructAdaptor += " + " + nameAdaptor;
        }
      }
      oldId = id;
    }

    //
    // Analyze counter group Construction
    //

    try {

      FileWriter writerStat = new FileWriter(this.nameOutputFastq.getParent()
          + "/differents_constructions_of_RT_adaptor_on_outliers.txt");

      writerStat.write(header + "\n\n");

      // get all adaptor construction in the localReporterNumberTimesAdaptor
      // object
      for (String constructionAdaptor : localReporterNumberTimesAdaptor
          .getCounterNames("Construction")) {
        writerStat.write(constructionAdaptor
            + " :  " + localReporterNumberTimesAdaptor
                .getCounterValue("Construction", constructionAdaptor)
            + "\n");
      }

      writerStat.close();

    } catch (Exception e) {

      e.printStackTrace();
    }

  }

  //
  // Main execution
  //

  /**
   * Method of the class CutadaptTrimmer to trim with the Trimmer interface.
   * @param leftLengthOutlier , the length of the left outlier
   * @param rightLengthOutlier , the length of the left outlier
   * @param sequence, the sequence of the read
   * @param id , the id of the read
   * @param quality , the quality of the read
   */
  public void preProcessSequence(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality) throws IOException {

    // write outliers in a fasta file
    UtilsTrimming.writeOutliers(leftLengthOutlier, rightLengthOutlier, sequence,
        id, new FastaWriter(fastaLeftOutlierFile),
        new FastaWriter(fastaRightOutlierFile));

  }

  /**
   * Method of the class CutadaptTrimmer to trim with the Trimmer interface
   */
  public void trim() {

    System.out.println("Begin use cutadapt !");

    String strandLeft = "-g";
    String strandRight = "-a";

    // Cutadapt execution for the left outlier
    cutadaptTrimming(this.fastaLeftOutlierFile, strandLeft,
        this.infoTrimLeftFile, this.outputTrimLeftFastaFile);

    // Cutadapt execution for the right outlier
    cutadaptTrimming(this.fastaRightOutlierFile, strandRight,
        this.infoTrimRightFile, this.outputTrimRightFastaFile);

    // Merge the output form cutadapt
    mergeTrimOutlier();

  }
}
