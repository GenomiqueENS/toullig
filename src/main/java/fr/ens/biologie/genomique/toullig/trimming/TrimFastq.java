package fr.ens.biologie.genomique.toullig.trimming;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder.OutlierPositionFinder;
import fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder.OutlierPositionFinderFactory;
import fr.ens.biologie.genomique.toullig.trimming.Trimmer.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

/**
 * Class to execute the cutadaptTrimming of ONT read. Created by birer on
 * 24/02/17.
 * @author Aurelien Birer
 */
public class TrimFastq {

  private String adaptorRT;
  private String adaptorStrandSwitching;

  private File samFile;
  private File fastqFile;
  private File workDir;
  private File outputFastqFile;
  private File adaptorFile;

  private String trimmer;
  private String mode;

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
   * @param outputFastqFile, a output fasqt directory
   * @param workDir, a work directory
   * @param trimmer, a trimmer to use
   * @param mode, a outlierpositionfinder to use
   * @throws IOException if an IO error occur
   */
  public TrimFastq(File samFile, File fastqFile, File adaptorFile,
      File outputFastqFile, File workDir, String trimmer, String mode)
      throws IOException {

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
    if (!outputFastqFile.exists()) {
      try {
        PrintWriter nameOutputFastqWriter =
            new PrintWriter(outputFastqFile.toString(), "UTF-8");
        nameOutputFastqWriter.close();
      } catch (IOException e) {
        e.printStackTrace();
      }

      this.outputFastqFile = outputFastqFile;
    } else {
      this.outputFastqFile = outputFastqFile;
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

    // test if trimmer is null
    if (trimmer != null) {

      this.trimmer = trimmer;

    } else {

      System.out.println("The trimmer is null !");
      System.exit(0);

    }

    // test if mode is null
    if (mode != null) {

      this.mode = mode;

    } else {

      System.out.println("The mode is null !");
      System.exit(0);

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
  private static void readSamFile(InputStream samInputStream,
      Map<String, InformationRead> workTrimmingMap) {

    // opne the sam File
    try (final SamReader inputSam = SamReaderFactory.makeDefault()
        .open(SamInputResource.of(samInputStream))) {

      int lengthBeginOutlier = 0;
      int lengthEndOutlier = 0;
      String quality = "";
      String sequence = "";

      // create a hash for multi-mapped read
      Set<String> multiMappedReadsSet = null;

      try {

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
          if (workTrimmingMap.containsKey(id)) {

            // add this read to the multi-mapped reads
            multiMappedReadsSet.add(id);

            // get the informationRead of existing read on the work trim map
            InformationRead informationRead = workTrimmingMap.get(id);

            // Select the largest alignement for a multi-mapped read
            if (informationRead.cigarLength <= cigarLength) {

              // add the read to the work Trimming map
              workTrimmingMap.put(id,
                  new InformationRead(sequence, quality, cigar,
                      lengthBeginOutlier, lengthEndOutlier, qFlag,
                      cigarLength));
            }

          } else {

            // add the read to the work Trimming map
            workTrimmingMap.put(id,
                new InformationRead(sequence, quality, cigar,
                    lengthBeginOutlier, lengthEndOutlier, qFlag, cigarLength));
          }

        }

      } catch (Exception e) {
        e.printStackTrace();
      }

      if (multiMappedReadsSet != null) {

        getLogger().info(
            "Number of multi-mapped reads: " + multiMappedReadsSet.size());

      }

      inputSam.close();

    } catch (IOException e) {

      e.printStackTrace();
    }

  }

  //
  // Set
  //

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

  public void setProcessStatsCutadapt(boolean processStatsCutadapt) {
    this.processStatsCutadapt = processStatsCutadapt;
  }

  //
  // Main execution
  //

  /**
   * Method of the class TrimFastq to execute the cutadaptTrimming.
   * @throws IOException if an IO error occur
   */
  public void execution() throws IOException {

    Map<String, InformationRead> workTrimmingMap = new HashMap<>();

    getLogger().info("add_index: " + this.addIndexOutlier);

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
    readSamFile(samInputStream, workTrimmingMap);

    // Declare the left outlier fasta for cutadapt
    File fastaLeftOutlierFile =
        new File(this.workDir + "/fastaFileLeftOutlier.fasta");

    // Declare the right outlier fasta for cutadapt
    File fastaRightOutlierFile =
        new File(this.workDir + "/fastaFileRightOutlier.fasta");

    // Declare the left info on trim for cutadapt
    File infoTrimLeftFile =
        new File(this.workDir + "/logCutadaptLeftOutlier.txt");

    // Declare the right info on trim for cutadapt
    File infoTrimRightFile =
        new File(this.workDir + "/logCutadaptRightOutlier.txt");

    // Create the OutlierPositionFinder Object with the correct method to the
    // OutlierPositionFactory
    OutlierPositionFinder outlierPositionFinder =
        OutlierPositionFinderFactory.newOutlierPositionFinder(workTrimmingMap,
            this.addIndexOutlier, this.fastqFile, this.lengthWindowSideWindow,
            this.thresholdSideWindow, this.mode);

    // Create the Trimmer Object with the correct method to the TrimmerFactory
    Trimmer trimmer = TrimmerFactory.newTrimmer(workTrimmingMap,
        this.outputFastqFile, outputTrimLeftFastaFile, outputTrimRightFastaFile,
        this.adaptorRT, this.adaptorStrandSwitching, this.errorRateCutadapt,
        fastaLeftOutlierFile, fastaRightOutlierFile, infoTrimLeftFile,
        infoTrimRightFile, this.adaptorFile, this.seedMismatchesTrimmomatic,
        this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold,
        this.trimmer, this.processStatsCutadapt);

    // execute the outlier position finder
    outlierPositionFinder.findOutliers(fastaLeftOutlierFile,
        fastaRightOutlierFile, trimmer);

    // execute the trimming
    if (trimmer != null) {
      trimmer.trim();
    }

  }
}