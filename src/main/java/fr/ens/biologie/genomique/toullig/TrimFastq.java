package fr.ens.biologie.genomique.toullig;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.toullig.trimming.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Class to execute the cutadaptTrimming of ONT read. Created by birer on
 * 24/02/17.
 */
public class TrimFastq implements AutoCloseable {

  private final Map<String, String[]> workTrimmingMap = new HashMap<>();
  private String adaptorRT;
  private String adaptorStrandSwitching;

  private File outputTrimLeftFastaFile;
  private File outputTrimRightFastaFile;

  private File samFile;
  private File fastqFile;
  private File workDir;
  private File nameOutputFastq;
  private File adaptorFile;

  private final boolean processPerfectTrim = true;
  private boolean processSideWindowTrim = false;
  private boolean processCutadapt = true;
  private boolean processTrimmomatic = false;
  private boolean processStatsCutadapt = false;

  private int addIndexOutlier = 15;
  private int lengthWindowSideWindow = 15;
  private double thresholdSideWindow = 0.8;
  private double errorRateCutadapt = 0.5;
  private int seedMismatchesTrimmomatic = 17;
  private int palindromeClipThresholdTrimmomatic = 30;
  private int simpleClipThreshold = 7;

  /**
   * Constructor of the class TrimFastq.
   * @param samFile, a sam File
   * @param fastqFile, a fastq File
   * @param adaptorFile, a adaptor (.txt) File
   * @param nameOutputFastq, a output fasqt directory
   * @param workDir, a work directory
   * @throws IOException if an IO error occur
   */
  public TrimFastq(File samFile, File fastqFile, File adaptorFile,
      File nameOutputFastq, File workDir) throws IOException {

    if (!samFile.exists()) {
      throw new IOException("The file " + samFile + " dont exist!");
    } else {
      this.samFile = samFile;
    }

    if (!fastqFile.exists()) {
      throw new IOException("The file " + fastqFile + " dont exist!");
    } else {
      this.fastqFile = fastqFile;
    }

    if (!nameOutputFastq.exists()) {
      try {
        PrintWriter nameOutputFastqWriter =
            new PrintWriter(nameOutputFastq.toString(), "UTF-8");
        nameOutputFastqWriter.close();
      } catch (IOException e) {
        e.printStackTrace();
      }

      this.nameOutputFastq = nameOutputFastq;
    } else {
      this.nameOutputFastq = nameOutputFastq;
    }

    if (!adaptorFile.exists()) {
      throw new IOException("The file " + adaptorFile + " dont exist!");
    } else {
      this.adaptorFile = adaptorFile;
    }

    if (!workDir.exists()) {
      throw new IOException("The directory " + workDir + " dont exist!");
    } else {
      this.workDir = workDir;
    }
  }

  //
  // Read Files
  //

  /**
   * Method of the class TrimFastq to read and get the adaptor in the adaptor
   * file.
   * @param adaptorFile, the file who contains adaptor
   */
  private void readAdaptorRTFile(File adaptorFile) {

    try (BufferedReader adaptorBufferedReader =
        new BufferedReader(new FileReader(adaptorFile))) {
      int i = 0;
      String line;
      while ((line = adaptorBufferedReader.readLine()) != null) {
        if (i == 1) {
          this.adaptorRT = line;
        }
        if (i == 3) {
          this.adaptorStrandSwitching = line;
        }
        i++;
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Method of the class TrimFastq to read and get the ID,sequence, quality and
   * CIGAR of a sam file.
   * @param samInputStream, the input stream file sam
   */
  private void readSamFile(InputStream samInputStream) {

    try (final SamReader inputSam = SamReaderFactory.makeDefault()
        .open(SamInputResource.of(samInputStream))) {

      String lengthBeginOutlier = "";
      String lengthEndOutlier = "";
      String quality = "";
      String sequence = "";

      HashSet<String> fastqHashMultimapped = new HashSet<>();

      for (SAMRecord samRecord : inputSam) {

        String cigar = samRecord.getCigarString();

        int qFlag = samRecord.getFlags();
        String id = samRecord.getReadName();
        int cigarLength = samRecord.getCigarLength();

        if (this.workTrimmingMap.containsKey(id)) {
          fastqHashMultimapped.add(id);
          String[] tabvalue = this.workTrimmingMap.get(id);
          if (Integer.parseInt(tabvalue[6]) <= cigarLength) {
            this.workTrimmingMap.put(id,
                new String[] {sequence, quality, cigar, lengthBeginOutlier,
                    lengthEndOutlier, "" + qFlag, "" + cigarLength});
          }

        } else {
          this.workTrimmingMap.put(id,
              new String[] {sequence, quality, cigar, lengthBeginOutlier,
                  lengthEndOutlier, "" + qFlag, "" + cigarLength});
        }
      }
      System.out.println(
          "Number of multi-mapped reads: " + fastqHashMultimapped.size());
      inputSam.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Method of the class TrimFastq to read a fastq file.
   * @param fastqFile, a BufferedReader fastq file
   */
  private void readFastqFile(File fastqFile) {

    try (FastqReader reader = new FastqReader(fastqFile)) {
      for (ReadSequence read : reader) {

        String header = read.getName();
        String[] part = header.split(" ");
        String ID = part[0];
        String sequence = read.getSequence();
        String quality = read.getQuality();
        String[] tabValue = this.workTrimmingMap.get(ID);
        tabValue[0] = sequence;
        tabValue[1] = quality;
      }
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  //
  // Set
  //

  /**
   * Method of the class TrimFastq to set the process SW in trim mode.
   */
  public void setProcessSideWindowTrim() {
    this.processSideWindowTrim = true;
    this.processCutadapt = false;
  }

  /**
   * Method of the class TrimFastq to set the process trimmomatic.
   */
  public void setProcessTrimmomatic() {
    this.processTrimmomatic = true;
    this.processCutadapt = false;
  }

  /**
   * Method of the class TrimFastq to set the process stats.
   */
  public void setProcessStats() {
    this.processStatsCutadapt = true;
  }

  /**
   * Method of the class TrimFastq to set the number of bases to add for the
   * outliers.
   * @param addIndexOutlier, a int
   */
  public void setAddIndexOutlier(int addIndexOutlier) {
    this.addIndexOutlier = addIndexOutlier;
  }

  /**
   * Method of the class TrimFastq to set the the length of the SW window.
   * @param lengthWindowSideWindow, a int
   */
  public void setLengthWindowSideWindow(int lengthWindowSideWindow) {
    this.lengthWindowSideWindow = lengthWindowSideWindow;
  }

  /**
   * Method of the class TrimFastq to set the treshold of the SW mode.
   * @param thresholdSideWindow, a double
   */
  public void setThresholdSideWindow(double thresholdSideWindow) {
    this.thresholdSideWindow = thresholdSideWindow;
  }

  /**
   * Method of the class TrimFastq to set the error rate for cutadapt.
   * @param errorRateCutadapt, a double
   */
  public void setErrorRateCutadapt(double errorRateCutadapt) {
    this.errorRateCutadapt = errorRateCutadapt;
  }

  /**
   * Method of the class TrimFastq to set the seed mismatches for trimmomatic.
   * @param seedMismatchesTrimmomatic, a int
   */
  public void setSeedMismatchesTrimmomatic(int seedMismatchesTrimmomatic) {
    this.seedMismatchesTrimmomatic = seedMismatchesTrimmomatic;
  }

  /**
   * Method of the class TrimFastq to set the palindrome clip treshold for
   * trimmomatic.
   * @param palindromeClipThresholdTrimmomatic, a int
   */
  public void setPalindromeClipThresholdTrimmomatic(
      int palindromeClipThresholdTrimmomatic) {
    this.palindromeClipThresholdTrimmomatic =
        palindromeClipThresholdTrimmomatic;
  }

  /**
   * Method of the class TrimFastq to set the simple clip treshold for
   * trimmomatic.
   * @param simpleClipThreshold, a int
   */
  public void setSimpleClipThreshold(int simpleClipThreshold) {
    this.simpleClipThreshold = simpleClipThreshold;
  }

  //
  // Main execution
  //

  /**
   * Method of the class TrimFastq to execute the cutadaptTrimming.
   * @throws IOException if an IO error occur
   * @throws InterruptedException if Interrupted IO error occur
   */
  public void execution() throws IOException {

    this.outputTrimLeftFastaFile =
        new File(this.workDir + "/outputFastaFileLeftOutlier.fastq");
    this.outputTrimRightFastaFile =
        new File(this.workDir + "/outputFastaFileRightOutlier.fastq");

    // Problem with ONT skip read for the RT adaptor (to many TTTT..)
    readAdaptorRTFile(this.adaptorFile);

    InputStream samInputStream = new FileInputStream(this.samFile);
    readSamFile(samInputStream);

    readFastqFile(this.fastqFile);

    OutlierPositionFinder outlierPositionFinder;
    Trimmer trimmer;

    File fastaLeftOutlierFile =
        new File(this.workDir + "/fastaFileLeftOutlier.fasta");
    File fastaRightOutlierFile =
        new File(this.workDir + "/fastaFileRightOutlier.fasta");

    File infoTrimLeftFile =
        new File(this.workDir + "/logCutadaptLeftOutlier.txt");
    File infoTrimRightFile =
        new File(this.workDir + "/logCutadaptRightOutlier.txt");

    if (this.processPerfectTrim) {
      outlierPositionFinder = new PerfectOutlierPositionFinder(
          this.workTrimmingMap, this.addIndexOutlier);
    } else {
      outlierPositionFinder =
          new SideWindowOutlierPositionFinder(this.lengthWindowSideWindow,
              this.thresholdSideWindow, this.workTrimmingMap);
    }

    if (this.processCutadapt) {
      trimmer = new CutadaptTrimmer(this.workTrimmingMap, this.nameOutputFastq,
          this.outputTrimLeftFastaFile, this.outputTrimRightFastaFile,
          this.adaptorRT, this.adaptorStrandSwitching, this.errorRateCutadapt,
          fastaLeftOutlierFile, fastaRightOutlierFile, infoTrimLeftFile,
          infoTrimRightFile);

    } else {
      trimmer = new TrimmomaticTrimmer(this.adaptorFile, this.nameOutputFastq,
          this.seedMismatchesTrimmomatic,
          this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold);
    }

    outlierPositionFinder.findOutliers(fastaLeftOutlierFile,
        fastaRightOutlierFile, trimmer);

    trimmer.trimming();

    if (processStatsCutadapt) {
      CutadaptTrimmer trimmingCutadapt =
          new CutadaptTrimmer(this.workTrimmingMap, this.nameOutputFastq,
              this.outputTrimLeftFastaFile, this.outputTrimRightFastaFile,
              this.adaptorRT, this.adaptorStrandSwitching,
              this.errorRateCutadapt, fastaLeftOutlierFile,
              fastaRightOutlierFile, infoTrimLeftFile, infoTrimRightFile);
      System.out.println("Start stat Left outlier");
      trimmingCutadapt.statsLogCutadapt(infoTrimLeftFile);
      System.out.println("Start stat Right outlier");
      trimmingCutadapt.statsLogCutadapt(infoTrimRightFile);
    }
  }

  @Override
  public void close() throws Exception {

  }
}