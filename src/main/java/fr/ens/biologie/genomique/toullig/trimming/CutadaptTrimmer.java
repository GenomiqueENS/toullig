package fr.ens.biologie.genomique.toullig.trimming;

import com.google.common.base.Strings;
import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.Sequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import fr.ens.biologie.genomique.toullig.Utils;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import static fr.ens.biologie.genomique.eoulsan.bio.Alphabets.AMBIGUOUS_DNA_ALPHABET;
import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;
import static fr.ens.biologie.genomique.toullig.Utils.complement;

/**
 * Class to execute the Cutadapt trimmer. Created by birer on 27/03/17.
 */
public class CutadaptTrimmer {

  private Map<String, String[]> workTrimmingMap;
  private File nameOutputFastq;
  private File outputTrimLeftFastaFile;
  private File outputTrimRightFastaFile;
  private String adaptorRetroTranscritpion;
  private String adaptorStrandSwitching;
  private Alphabet alphabet = AMBIGUOUS_DNA_ALPHABET;
  private double errorRateCutadapt;

  public CutadaptTrimmer(Map<String, String[]> workTrimmingMap,
      File nameOutputFastq, File outputTrimLeftFastaFile,
      File outputTrimRightFastaFile, String adaptorRetroTranscritpion,
      String adaptorStrandSwitching, double errorRateCutadapt) {

    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = nameOutputFastq;
    this.outputTrimLeftFastaFile = outputTrimLeftFastaFile;
    this.outputTrimRightFastaFile = outputTrimRightFastaFile;
    this.adaptorRetroTranscritpion = adaptorRetroTranscritpion;
    this.adaptorStrandSwitching = adaptorStrandSwitching;
    this.errorRateCutadapt = errorRateCutadapt;
  }

  //
  // Merge outlier
  //

  /**
   * Method of the class CutadaptTrimmer to merge the results of cutadapt with
   * the trim sequence.
   * @throws IOException
   */
  public void mergeTrimOutlier() throws IOException {

    HashMap<String, String> fastaLeftHash = new HashMap<String, String>();
    HashMap<String, String> fastaRightHash = new HashMap<String, String>();

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
        int lengthBeginOutlier;
        int lengthEndOutlier;
        String[] tabValue = this.workTrimmingMap.get(id);
        String cigar = tabValue[2];
        String score = tabValue[1];
        String sequence = tabValue[0];
        String beginOutlier = tabValue[3];
        String endOutlier = tabValue[4];

        if (!cigar.equals("*")) {

          if (beginOutlier.equals("")) {
            lengthBeginOutlier = 0;
          } else {
            lengthBeginOutlier = Integer.parseInt(beginOutlier);
          }
          if (endOutlier.equals("")) {
            lengthEndOutlier = 0;
          } else {
            lengthEndOutlier = Integer.parseInt(endOutlier);
          }

          if (lengthBeginOutlier < 0) {
            lengthBeginOutlier = 0;
          }
          if (lengthEndOutlier < 0) {
            lengthEndOutlier = 0;
          }

          String mainSequenceWithoutOutlier = "";

          if ((lengthBeginOutlier + lengthEndOutlier) < sequence.length()) {
            mainSequenceWithoutOutlier = sequence.substring(lengthBeginOutlier,
                sequence.length() - lengthEndOutlier);
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

          int leftLength = lengthBeginOutlier - leftSequence.length();
          int rightLength =
              (sequence.length() - lengthEndOutlier) + rightSequence.length();

          if (lengthBeginOutlier - leftSequence.length() < 0) {
            leftLength = 0;
          }

          if (rightLength > score.length()) {
            rightLength = score.length();
          }

          // case with the addition of index (>0)
          if (leftLength > rightLength) {
            leftLength = rightLength;
          }

          String sequenceTrim = sequence.substring(leftLength, rightLength);

          String scoreTrim = score.substring(leftLength, rightLength);

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
   * @throws IOException
   * @throws InterruptedException
   */
  public void cutadaptTrim(File fastaOutlierFile, String strand,
      File infoTrimFile, File pathOutputTrimFasta)
      throws IOException, InterruptedException {

    try {
      Utils utils = new Utils();
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

      StringBuffer reverseAdaptorRTStringBuffer =
          new StringBuffer(this.adaptorRetroTranscritpion);
      StringBuffer reverseAdaptorSwithStrandStringBuffer =
          new StringBuffer(this.adaptorStrandSwitching);
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

      int exitproc = proc.waitFor();
      // if(!proc.waitFor(1, TimeUnit.MINUTES)) {
      // //timeout - kill the process.
      // proc.destroy(); // consider using destroyForcibly instead
      // }
    } catch (IOException e) {
      e.printStackTrace(); // or log it, or otherwise handle it
    } catch (InterruptedException ie) {
      ie.printStackTrace(); // or log it, or otherwise handle it
    }
  }

  /**
   * Method of the class CutadaptTrimmer to display log of cutadapt (delete the
   * option --quiet).
   * @param proc, a Processus
   * @throws IOException
   */
  private void getLogCutadapt(Process proc) throws IOException {

    BufferedReader stdInput =
        new BufferedReader(new InputStreamReader(proc.getInputStream()));
    BufferedReader stdError =
        new BufferedReader(new InputStreamReader(proc.getErrorStream()));

    // read the output from the command
    System.out.println("Here is the standard output of the command:\n");
    String s = null;
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
   * Method of the class CutadaptTrimmer to make some stats of the trimming
   * @param infoTrimFile, a path to store stats in a file
   * @throws IOException
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
}
