package fr.ens.biologie.genomique.toullig;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.toullig.trimming.CutadaptTrimmer;
import fr.ens.biologie.genomique.toullig.trimming.PerfectOutlierPositionFinder;
import fr.ens.biologie.genomique.toullig.trimming.SideWindowOutlierPositionFinder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Class to execute the trimming of ONT read. Created by birer on 24/02/17.
 */
public class TrimFastq implements AutoCloseable {

  private Map<String, String[]> workTrimmingMap =
      new HashMap<String, String[]>();
  private String adaptorRT;
  private String adaptorStrandSwitching;

  private File outputTrimLeftFastaFile;
  private File outputTrimRightFastaFile;

  private File samFile;
  private File fastqFile;
  private File workDir;
  private File nameOutputFastq;
  private File adaptorFile;

  private boolean processPerfectTrim = true;
  private boolean processSideWindowTrim = false;
  private boolean processCutadapt = true;
  private boolean processTrimmomatic = false;
  private boolean processStats = false;

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
   * @throws IOException
   * @throws InterruptedException
   */
  public TrimFastq(File samFile, File fastqFile, File adaptorFile,
      File nameOutputFastq, File workDir)
      throws IOException, InterruptedException {

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
   * @throws IOException
   */
  private void readAdaptorRTFile(File adaptorFile) throws IOException {

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
   * @throws IOException
   */
  private void readSamFile(InputStream samInputStream) throws IOException {

    try (final SamReader inputSam = SamReaderFactory.makeDefault()
        .open(SamInputResource.of(samInputStream))) {

      String lengthBeginOutlier = "";
      String lengthEndOutlier = "";
      String quality = "";
      String sequence = "";

      HashSet<String> fastqHashMultimapped = new HashSet<String>();

      for (SAMRecord samRecord : inputSam) {

        String cigar = samRecord.getCigarString();

        int qFlag = samRecord.getFlags();
        String ID = samRecord.getReadName();
        int cigarLength = samRecord.getCigarLength();

        if (this.workTrimmingMap.containsKey(ID)) {
          fastqHashMultimapped.add(ID);
          String[] tabvalue = this.workTrimmingMap.get(ID);
          if (Integer.parseInt(tabvalue[6]) <= cigarLength) {
            this.workTrimmingMap.put(ID,
                new String[] {sequence, quality, cigar, lengthBeginOutlier,
                    lengthEndOutlier, "" + qFlag, "" + cigarLength});
          } else {
            continue;
          }

        } else {
          this.workTrimmingMap.put(ID,
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
   * @throws IOException
   */
  private void readFastqFile(File fastqFile) throws IOException {

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
   * Method of the class TrimFastq to set the process P in trim mode.
   * @param processPerfectTrim, a boolean
   */
  public void setProcessPerfectTrim(boolean processPerfectTrim) {
    this.processPerfectTrim = processPerfectTrim;
  }

  /**
   * Method of the class TrimFastq to set the process SW in trim mode.
   */
  public void setProcessSideWindowTrim() {
    this.processSideWindowTrim = true;
  }

  /**
   * Method of the class TrimFastq to set the process cutadapt.
   * @param processCutadapt, a boolean
   */
  public void setProcessCutadapt(boolean processCutadapt) {
    this.processCutadapt = processCutadapt;
  }

  /**
   * Method of the class TrimFastq to set the process trimmomatic.
   */
  public void setProcessTrimmomatic() {
    this.processTrimmomatic = true;
  }

  /**
   * Method of the class TrimFastq to set the process stats.
   */
  public void setProcessStats() {
    this.processStats = true;
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
  // Execution
  //

  /**
   * Method of the class TrimFastq to execute the process with cutadapt.
   * @param fastqFile, a fastq file
   * @throws IOException
   * @throws InterruptedException
   */
  private void executionTrimWithCutadapt(File fastqFile)
      throws IOException, InterruptedException {

    File fastaLeftOutlierFile =
        new File(this.workDir + "/fastaFileLeftOutlier.fasta");
    File fastaRightOutlierFile =
        new File(this.workDir + "/fastaFileRightOutlier.fasta");

    File infoTrimLeftFile =
        new File(this.workDir + "/logCutadaptLeftOutlier.txt");
    File infoTrimRightFile =
        new File(this.workDir + "/logCutadaptRightOutlier.txt");

    readFastqFile(fastqFile);

    System.out.println("Trim begin !");

    //

    if (this.processPerfectTrim) {
      PerfectOutlierPositionFinder perfectOutlierPositionFinder =
          new PerfectOutlierPositionFinder(this.workTrimmingMap,
              this.nameOutputFastq, this.addIndexOutlier, this.adaptorFile,
              this.seedMismatchesTrimmomatic,
              this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold,
              this.processCutadapt, this.processTrimmomatic);

      perfectOutlierPositionFinder.findOutliers(fastaLeftOutlierFile,
          fastaRightOutlierFile);
    }
    if (this.processSideWindowTrim) {
      SideWindowOutlierPositionFinder sideWindowOutlierPositionFinder =
          new SideWindowOutlierPositionFinder(this.lengthWindowSideWindow,
              this.thresholdSideWindow, this.workTrimmingMap,
              this.nameOutputFastq, this.adaptorFile,
              this.seedMismatchesTrimmomatic,
              this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold,
              this.processCutadapt, this.processTrimmomatic);

      sideWindowOutlierPositionFinder.getOutliers(fastaLeftOutlierFile,
          fastaRightOutlierFile);
    }

    CutadaptTrimmer trimmingCutadapt = new CutadaptTrimmer(this.workTrimmingMap,
        this.nameOutputFastq, this.outputTrimLeftFastaFile,
        this.outputTrimRightFastaFile, this.adaptorRT,
        this.adaptorStrandSwitching, this.errorRateCutadapt);

    System.out.println("Begin use cutadapt !");

    String strandLeft = "-g";
    String strandRight = "-a";

    // Cutadapt execution
    trimmingCutadapt.cutadaptTrim(fastaLeftOutlierFile, strandLeft,
        infoTrimLeftFile, this.outputTrimLeftFastaFile);
    trimmingCutadapt.cutadaptTrim(fastaRightOutlierFile, strandRight,
        infoTrimRightFile, this.outputTrimRightFastaFile);

    // Merge the output form cutadapt
    trimmingCutadapt.mergeTrimOutlier();

    if (processStats) {
      System.out.println("Start stat Left outlier");
      trimmingCutadapt.statsLogCutadapt(infoTrimLeftFile);
      System.out.println("Start stat Right outlier");
      trimmingCutadapt.statsLogCutadapt(infoTrimRightFile);
    }

  }

  /**
   * Method of the class TrimFastq to execute the process with trimmomatic.
   * @param fastqFile, a fastq file
   * @throws IOException
   * @throws InterruptedException
   */
  private void executionTrimWithTrimmomatic(File fastqFile)
      throws IOException, InterruptedException {

    readFastqFile(fastqFile);

    System.out.println("Begin use trimmomatic !");

    if (this.processPerfectTrim) {
      PerfectOutlierPositionFinder perfectOutlierPositionFinder =
          new PerfectOutlierPositionFinder(this.workTrimmingMap,
              this.nameOutputFastq, this.addIndexOutlier, this.adaptorFile,
              this.seedMismatchesTrimmomatic,
              this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold,
              this.processCutadapt, this.processTrimmomatic);

      perfectOutlierPositionFinder.findOutliers(new File(""), new File(""));
    }
    if (this.processSideWindowTrim) {
      SideWindowOutlierPositionFinder sideWindowOutlierPositionFinder =
          new SideWindowOutlierPositionFinder(this.lengthWindowSideWindow,
              this.thresholdSideWindow, this.workTrimmingMap,
              this.nameOutputFastq, this.adaptorFile,
              this.seedMismatchesTrimmomatic,
              this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold,
              this.processCutadapt, this.processTrimmomatic);

      sideWindowOutlierPositionFinder.getOutliers(new File(""), new File(""));
    }
  }

  //
  // Main execution
  //

  /**
   * Method of the class TrimFastq to execute the trimming.
   * @throws IOException
   * @throws InterruptedException
   */
  public void execution() throws IOException, InterruptedException {

    this.outputTrimLeftFastaFile =
        new File(this.workDir + "/outputFastaFileLeftOutlier.fastq");
    this.outputTrimRightFastaFile =
        new File(this.workDir + "/outputFastaFileRightOutlier.fastq");

    // Problem with ONT skip read for the RT adaptor (to many TTTT..)
    readAdaptorRTFile(this.adaptorFile);

    InputStream samInputStream = new FileInputStream(this.samFile);
    readSamFile(samInputStream);

    if (this.processCutadapt) {
      executionTrimWithCutadapt(this.fastqFile);
    }
    if (this.processTrimmomatic) {
      executionTrimWithTrimmomatic(this.fastqFile);
    }
  }

  @Override
  public void close() throws Exception {

  }
}
