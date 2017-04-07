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

  private final Map<String, InformationRead> workTrimmingMap = new HashMap<>();
  private String adaptorRT;
  private String adaptorStrandSwitching;

  private File samFile;
  private File fastqFile;
  private File workDir;
  private File nameOutputFastq;
  private File adaptorFile;

  private boolean processCutadapt = true;
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

    // test if the sam File exist
    if (!samFile.exists()) {
      throw new IOException("The file " + samFile + " dont exist!");
    } else {
      this.samFile = samFile;
    }

    // test if the fastq File exist
    if (!fastqFile.exists()) {
      throw new IOException("The file " + fastqFile + " dont exist!");
    } else {
      this.fastqFile = fastqFile;
    }

    // test if the fastq output File exist
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

    // test if the adaptor File exist
    if (!adaptorFile.exists()) {
      throw new IOException("The file " + adaptorFile + " dont exist!");
    } else {
      this.adaptorFile = adaptorFile;
    }

    // test if the working Directory exist
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

    // open the adaptor File
    try (BufferedReader adaptorBufferedReader =
        new BufferedReader(new FileReader(adaptorFile))) {

      int i = 0;
      String line;

      // read the adaptor File
      while ((line = adaptorBufferedReader.readLine()) != null) {

        // get the first adaptor sequence
        if (i == 1) {
          this.adaptorRT = line;
        }

        // get the second adaptor sequence
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

    // opne the sam File
    try (final SamReader inputSam = SamReaderFactory.makeDefault()
        .open(SamInputResource.of(samInputStream))) {

      int lengthBeginOutlier = 0;
      int lengthEndOutlier = 0;
      String quality = "";
      String sequence = "";

      // create a hash for multi-mapped read
      HashSet<String> fastqHashMultimapped = new HashSet<>();

      // read sam File
      for (SAMRecord samRecord : inputSam) {

        // Get cigar information
        String cigar = samRecord.getCigarString();

        // Get qFlag information
        int qFlag = samRecord.getFlags();

        // Get id information
        String id = samRecord.getReadName();

        // Get cigarLength information
        int cigarLength = samRecord.getCigarLength();

        // test if the work Map already contains this id
        if (this.workTrimmingMap.containsKey(id)) {

          // add this read to the multi-mapped reads
          fastqHashMultimapped.add(id);

          // get the informationRead of existing read on the work trimming map
          InformationRead informationRead = this.workTrimmingMap.get(id);

          // Select the largest alignement for a multi-mapped read
          if (informationRead.cigarLength <= cigarLength) {

            // add the read to the work Trimming map
            this.workTrimmingMap.put(id,
                new InformationRead(sequence, quality, cigar,
                    lengthBeginOutlier, lengthEndOutlier, qFlag, cigarLength));
          }

        } else {

          // add the read to the work Trimming map
          this.workTrimmingMap.put(id, new InformationRead(sequence, quality,
              cigar, lengthBeginOutlier, lengthEndOutlier, qFlag, cigarLength));
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

    // open the fastq File
    try (FastqReader reader = new FastqReader(fastqFile)) {

      // read the fastq file
      for (ReadSequence read : reader) {

        // get header
        String header = read.getName();

        // split header for get id
        String[] part = header.split(" ");

        // get id
        String id = part[0];

        // get sequence
        String sequence = read.getSequence();

        // get quality
        String quality = read.getQuality();

        // get the information read of the id corresponding in the work trimming
        // map
        InformationRead informationRead = this.workTrimmingMap.get(id);

        // set the sequence of the read
        informationRead.sequence = sequence;

        // set the quality of the read
        informationRead.quality = quality;
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
    boolean processSideWindowTrim = true;
    this.processCutadapt = false;
  }

  /**
   * Method of the class TrimFastq to set the process trimmomatic.
   */
  public void setProcessTrimmomatic() {
    boolean processTrimmomatic = true;
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
   */
  public void execution() throws IOException {

    System.out.println("add_index: " + this.addIndexOutlier);
    System.out.println("error rate cutadapt: " + this.errorRateCutadapt);

    // Declare the left outlier output fasta for cutadapt
    File outputTrimLeftFastaFile =
        new File(this.workDir + "/outputFastaFileLeftOutlier.fastq");

    // Declare the right outlier output fasta for cutadapt
    File outputTrimRightFastaFile =
        new File(this.workDir + "/outputFastaFileRightOutlier.fastq");

    // Problem with ONT skip read for the RT adaptor (to many TTTT..)
    // read the adaptor File
    readAdaptorRTFile(this.adaptorFile);

    // create InputStream for read sam File
    InputStream samInputStream = new FileInputStream(this.samFile);

    // read the sam File
    readSamFile(samInputStream);

    // read the fastq File
    readFastqFile(this.fastqFile);

    // call the OutlierPositionFinder interface
    OutlierPositionFinder outlierPositionFinder;

    // call the Trimmer interface
    Trimmer trimmer;

    // Declare the left outlier fasta for cutadapt
    File fastaLeftOutlierFile =
        new File(this.workDir + "/fastaFileLeftOutlier.fasta");

    // Declare the right outlier fasta for cutadapt
    File fastaRightOutlierFile =
        new File(this.workDir + "/fastaFileRightOutlier.fasta");

    // Declare the left info on trimming for cutadapt
    File infoTrimLeftFile =
        new File(this.workDir + "/logCutadaptLeftOutlier.txt");

    // Declare the right info on trimming for cutadapt
    File infoTrimRightFile =
        new File(this.workDir + "/logCutadaptRightOutlier.txt");

    // Set the perfect method on true
    boolean processPerfectTrim = true;

    // test to process the perfect method
    if (processPerfectTrim) {

      // call PerfectOutlierPositionFinder constructor
      outlierPositionFinder = new PerfectOutlierPositionFinder(
          this.workTrimmingMap, this.addIndexOutlier);
    } else {

      // call SideWindowOutlierPositionFinder constructor
      outlierPositionFinder =
          new SideWindowOutlierPositionFinder(this.lengthWindowSideWindow,
              this.thresholdSideWindow, this.workTrimmingMap);
    }

    // test to process the cutadapt trimmer
    if (this.processCutadapt) {

      // call CutadaptTrimmer constructor
      trimmer = new CutadaptTrimmer(this.workTrimmingMap, this.nameOutputFastq,
          outputTrimLeftFastaFile, outputTrimRightFastaFile, this.adaptorRT,
          this.adaptorStrandSwitching, this.errorRateCutadapt,
          fastaLeftOutlierFile, fastaRightOutlierFile, infoTrimLeftFile,
          infoTrimRightFile);

    } else {

      // call TrimmomaticTrimmer constructor
      trimmer = new TrimmomaticTrimmer(this.adaptorFile, this.nameOutputFastq,
          this.seedMismatchesTrimmomatic,
          this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold);
    }

    // execute the outlier position finder
    outlierPositionFinder.findOutliers(fastaLeftOutlierFile,
        fastaRightOutlierFile, trimmer);

    // execute the trimming
    trimmer.trimming();

    // test to execute stats for cutadapt trimming
    if (processStatsCutadapt) {

      // call CutadaptTrimmer constructor
      CutadaptTrimmer trimmingCutadapt = new CutadaptTrimmer(
          this.workTrimmingMap, this.nameOutputFastq, outputTrimLeftFastaFile,
          outputTrimRightFastaFile, this.adaptorRT, this.adaptorStrandSwitching,
          this.errorRateCutadapt, fastaLeftOutlierFile, fastaRightOutlierFile,
          infoTrimLeftFile, infoTrimRightFile);

      System.out.println("Start stat Left outlier");
      // execute for the left outlier the stats of cutadapt
      trimmingCutadapt.statsLogCutadapt(infoTrimLeftFile);

      System.out.println("Start stat Right outlier");
      // execute for the right outlier the stats of cutadapt
      trimmingCutadapt.statsLogCutadapt(infoTrimRightFile);
    }
  }

  /**
   * Method to close File of the Autoclosable class.
   * @throws Exception , if an exception occur
   */
  @Override
  public void close() throws Exception {

  }
}