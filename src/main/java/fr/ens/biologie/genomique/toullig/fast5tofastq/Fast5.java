package fr.ens.biologie.genomique.toullig.fast5tofastq;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 * This class read a Fast5 of a minION (ONT) run. It's a HDF5 format file.
 * @author Aurelien Birer
 */
public class Fast5 implements AutoCloseable {

  private final Version version;
  private final Type type;
  private final Status status;
  private final ChemistryVersion chemistryVersion;
  private final Basecaller basecaller;
  private final IHDF5Reader reader;
  private final File fast5File;

  /**
   * Constructor of the Fast5 class.
   * @param fast5File constructor
   */
  public Fast5(File fast5File) throws ParseException {

    this.fast5File = fast5File;
    this.reader = readFast5File(fast5File);
    this.status = readStatus();
    this.basecaller = readBasecaller();
    this.version = readVersion();
    this.type = readType();
    this.chemistryVersion = readChemistryVersion();

  }

  //
  // Main declarations
  //

  /**
   * Close the fast5 file.
   */
  public void close() {
    this.reader.close();
  }

  /**
   * Values of the variable Version that design the version of the fast5 format.
   */
  public enum Version {
    V1_0, V1_1 // The version of the fast5 file ( 1.1 is here since 11/2015)
  }

  /**
   * Values of the variable Type that design the type of experimental design.
   */
  public enum Type {
    TYPE_1D, TYPE_2D, TYPE_1D2 // Type of the sequencing
  }

  /**
   * Values of the variable Status that design the state of the fast5 file (if
   * it was basecalled or not).
   */
  public enum Status {
    PRE_BASECALLING, AFTER_BASECALLING // The status of the fast5 file
  }

  /**
   * Values of the variable ChemistryVersion that design the chemical kit use.
   */
  public enum ChemistryVersion {
    R7_3, R9, R9_4, R9_5 // Chemical version of the kit
  }

  /**
   * Values of the variable Basecaller that design the basecaller use for the basecalling.
   */
  public enum Basecaller {
    METRICHOR, ALBACORE
  }

  //
  // read files
  //

  /**
   * Method who use the class IDFH5Reader to read the fast5 file.
   * @return a hdf5 file open
   */

  private static IHDF5Reader readFast5File(File fast5File) {

    // Get the object to read a hdf5 file
    IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
    return hdf5Fac.openForReading(fast5File);
  }

  /**
   * Method of the variable Version who obtain the version of the fast5 format.
   * @return a version with the version of the fast5 format
   */
  private Version readVersion() {
    if (!isBasecalled()) {
      return null;
    }
    if (reader.isGroup("/Analyses/Basecall_2D_000/BaseCalled_template")) {
      return Version.V1_0;
    }
    if (reader.isGroup("/Analyses/Basecall_1D_000/BaseCalled_template")) {
      return Version.V1_1;
    }
    return null;
  }

  /**
   * Method of the variable Type who obtain the type of experimental design.
   * @return a type with the type of sequencing done
   */
  private Type readType() {
    if (!isBasecalled()) {
      return null;
    }
    if (reader.isGroup("/Analyses/Basecall_2D_000")) {
      return Type.TYPE_2D;
    } else {
      return Type.TYPE_1D;
    }
  }

  /**
   * Method of the variable Status who obtain the state of the fast5 file (if it
   * was basecalled or not).
   * @return a status with the status of the fast5 file
   */
  private Status readStatus() {

    // test if the reader is open
    if (this.reader.getFile() == null) {
      throw new IllegalStateException("The file is closed");
    }

    // test if the fast5 file is basecalled
    if (reader.isGroup("/Analyses/Basecall_1D_000")) {
      return Status.AFTER_BASECALLING;
    } else {
      return Status.PRE_BASECALLING;
    }
  }

  /**
   * Method of the variable ChemistryVersion who obtain the chemical kit use.
   * @return a ChemistryVersion with the chemical kit version use
   */

  private ChemistryVersion readChemistryVersion() {

    //
    // Case of not basecalled fast5 file
    //

    // Case of not basecalled fast5 file and contains a special group
    if (!isBasecalled() && reader.isGroup("/Raw/Reads")) {
      return ChemistryVersion.R9;
    }

    // Case of not basecalled fast5 file and contains a special group
    if (!isBasecalled()
            && reader.isGroup("/Analyses/EventDetection_000/Reads")) {
      return ChemistryVersion.R7_3;
    }

    //
    // Case of basecalled fast5 file
    //

    // test if the chemistry version is R7.3
    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r7.3_")) {
      return ChemistryVersion.R7_3;
    }

    // test if the chemistry version is R9
    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r9_")) {
      return ChemistryVersion.R9;

    }

    // test if the chemistry version is R9.4
    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r94_")) {
      return ChemistryVersion.R9_4;

    }
    return null;
  }


  /**
   * Method of the variable Basecaller who obtain the basecaller use.
   * @return a Basecaller with the basecaller use
   */
  private Basecaller readBasecaller(){

      //test if the file is basecalled
      if (!isBasecalled()) {
          return null;
      }

    // test if the basecaller is Metrichor by a specific Metrichor field
    if (reader.getStringAttribute("/Analyses/Basecall_1D_000","name").equals("ONT Sequencing Workflow")){
      return Basecaller.METRICHOR;
    }else{
        // test if the basecaller is Albacore by a specific Albacore field
        if (reader.getStringAttribute("/Analyses/Basecall_1D_000",
                "name").equals("ONT Albacore Sequencing Software")){
            return Basecaller.ALBACORE;
        }
      }


    return null;
  }

  //
  //
  // Important getters
  //
  //

  /**
   * Getter of the name of tha fast5 file.
   * @return a string of name of fast5 file
   */
  public String getNameFast5File() {
    return this.fast5File.toString();
  }

  /**
   * Getter of the variable Status.
   * @return a status with the status of the fast5 file
   */
  public Status getStatus() {
    return this.status;
  }

  /**
   * Getter of the variable Version.
   * @return a version with the version of the fast5 format
   */
  public Version getVersion() {
    return this.version;
  }

  /**
   * Getter of the variable Type.
   * @return a type with the type of sequencing done
   */
  public Type getType() {
    return this.type;
  }

  /**
   * Getter of the variable RVersion.
   * @return a rversion with the chemical kit version use
   */
  public ChemistryVersion getChemistryVersion() {
    return this.chemistryVersion;
  }

  //
  //
  // macro
  //
  //

  /**
   * Boolean shortcut to know if the file is barcoded.
   * @return a boolean with the barcoded information
   */
  public boolean isBarcoded() {
    return reader.isGroup("/Analyses/Barcoding_000");
  }

  /**
   * Boolean shortcut to know if the file is basecalled.
   * @return a boolean with the basecalled information
   */
  public boolean isBasecalled() {
    return this.status == Status.AFTER_BASECALLING;
  }

  /**
   * Boolean shortcut to know if the experiement is 2D (or 1D in the opposite
   * case).
   * @return a boolean with the type 2D information
   */
  public boolean is2D() {
    return this.type == Type.TYPE_2D;
  }

  //
  // UniqueGlobalkey-tracking_id Group Information
  //

  /**
   * Method of the class Fast5 to obtain the serial number of the minION in the
   * fast5 file.
   * @return a string with the serial number of the minION
   */
  public String getNumMinION() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "device_id");
  }

  /**
   * Method of the class Fast5 to obtain the flowcell id in the fast5 file.
   * @return a string with the flowcell id
   */
  public String getFlowcellId() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "flow_cell_id");
  }

  /**
   * Method of the class Fast5 to obtain the MinKnow version in the fast5 file.
   * @return a string with the MinKnow version
   */
  public String getMinknowVersion() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id", "version");
  }

  /**
   * Method of the class Fast5 to obtain the date of the sequencing in the fast5
   * file.
   * @return a date of the sequencing
   */
  public Date getdDateExp() throws ParseException {

      // This function cann return error because the problem of teh date format is related with MinKNOW version and i dont know the version of MinKNOW where they change the format of the Date.

    // test if the Chemistry version is R9_4
    if(this.chemistryVersion==ChemistryVersion.R9_4 && this.chemistryVersion==ChemistryVersion.R9_5){
        SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MMM-ddTHH:mm:ssZ");
        String date = reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
                "exp_start_time");
        return formatter.parse(date);
    }else{
        String dateInt = reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
                "exp_start_time");
        return new Date(Long.parseLong(dateInt) * 1000);
    }
  }

  /**
   * Method of the class Fast5 to obtain the protocol id in the fast5 file.
   * @return a string with the protocol id
   */
  public String getProtocolRunId() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "protocol_run_id");
  }

  /**
   * Method of the class Fast5 to obtain the hostname (experimenter) in the
   * fast5 file.
   * @return a string with the host name
   */
  public String getHostname() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "hostname");
  }

  /**
   * Method of the class Fast5 to obtain the Operating System in the fast5 file.
   * @return a string with the Operating System
   */
  public String getOS() {
    return reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "operating_system");
  }

  //
  // UniqueGlobalkey-context_tags Group Information
  //

  /**
   * Method of the class Fast5 to obtain the name of the experiment kit in the
   * fast5 file.
   * @return a string with the experiment kit
   */
  public String getExperimentKit() {
    return reader.getStringAttribute("/UniqueGlobalKey/context_tags",
            "experiment_kit");
  }

  /**
   * Method of the class Fast5 to obtain the name of the experiment type in the
   * fast5 file.
   * @return a string with the experiment type
   */
  public String getExperimentType() {
    if (reader.hasAttribute("/UniqueGlobalKey/context_tags",
            "experiment_type")) {
      return reader.getStringAttribute("/UniqueGlobalKey/context_tags",
              "experiment_type");
    }
    return null;
  }

  /**
   * Method of the class Fast5 to obtain the frequency of the sample in the
   * fast5 file.
   * @return an int with sample frequency
   */
  public int getSampleFrequency() {
    return Integer.parseInt(reader.getStringAttribute(
            "/UniqueGlobalKey/context_tags", "sample_frequency"));
  }

  /**
   * Method of the class Fast5 to obtain the channel number of the pore in the
   * fast5 file.
   * @return an int with sample frequency
   */
  public int getChannelNumber() {
    return Integer.parseInt(reader
            .getStringAttribute("/UniqueGlobalKey/channel_id", "channel_number"));
  }

  //
  //
  // Basecalling information getters
  //
  //

  /**
   * Method of the class Fast5 to obtain the number of the read in the fast5
   * file.
   * @return an int with number of the read
   */
  public int getNumberRead() {

    // test if the fast5 file is basecalled and is R9
    if (!isBasecalled() && getChemistryVersion() == ChemistryVersion.R9) {
      String reads = reader.getAllGroupMembers("/Raw/Reads").get(0);
      return Integer.parseInt(reads.substring(reads.indexOf('_') + 1));
    }

    // test if the fast5 file is basecalled and is R7.3
    if (!isBasecalled() && getChemistryVersion() == ChemistryVersion.R7_3) {
      String reads = reader
              .getAllGroupMembers("/Analyses/EventDetection_000/Reads").get(0);
      return Integer.parseInt(reads.substring(reads.indexOf('_') + 1));
    }
    return Integer.parseInt(reader.getStringAttribute(
            "/Analyses/Basecall_1D_000/Configuration/general", "read_id"));
  }

  /**
   * Method of the class Fast5 to obtain the version of the sub-modules of Metrichor
   * @return a string with the version
   */
  public String getSubModuleMetrichorVersion() {

      if(this.basecaller==Basecaller.METRICHOR){

          // test if the fast5 file is basecalled
          if (!isBasecalled()) {
              return null;
          }
          // get the version of chimaera
          String chimaeraVersion = reader
                  .getStringAttribute("/Analyses/Basecall_1D_000", "chimaera version");

          // get the version of dragonet
          String dragonetVersion = reader
                  .getStringAttribute("/Analyses/Basecall_1D_000", "dragonet version");
          return "chimaera v" + chimaeraVersion + " | dragonet v" + dragonetVersion;
      }

      return null;

  }

    /**
     * Method of the class Fast5 to obtain the version of Albacore
     * @return a string with the version of Albacore
     */
    public String getAlbacoreVersion() {

        if(this.basecaller==Basecaller.ALBACORE){

            // test if the fast5 file is basecalled
            if (!isBasecalled()) {
                return null;
            }

            // get the version of chimaera
            return reader
                    .getStringAttribute("/Analyses/Basecall_1D_000", "version");
        }

        return null;

    }

  /**
   * Method of the class Fast5 to obtain alignement of the complement and the
   * template sequence in the fast5 file.
   * @return a string with the alignement
   */
  public String getAlignment() {

    // test if the fast5 file is basecalled and is 2D
    if (!isBasecalled() || !is2D()) {
      return null;
    }
    return reader
            .readString("/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment");
  }

  /**
   * Method of the class Fast5 to obtain the length of the template sequence in
   * the fast5 file.
   * @return an int with the length of the template strand
   */
  public int getTemplateLength() {

    // test if the fast5 file is basecalled
    if (!isBasecalled()) {
      return 0;
    }

    return reader.getIntAttribute(
            "/Analyses/Basecall_1D_000/Summary/basecall_1d_template",
            "sequence_length");
  }

  /**
   * Method of the class Fast5 to obtain the length of the complemente sequence
   * in the fast5 file.
   * @return an int with the length of the complemente strand
   */
  public int getComplementeLength() {

    // test if the fast5 file is basecalled and is 2D
    if (!isBasecalled() || !is2D()) {
      return 0;
    }

    return reader.getIntAttribute(
            "/Analyses/Basecall_1D_000/Summary/basecall_1d_complement",
            "sequence_length");
  }

  /**
   * Method of the class Fast5 to obtain the serial number of the barcode in the
   * fast5 file.
   * @return a string with the barcode id
   */
  public String getNumBarcode() {

    // test if the fast5 file is basecalled
    if (!isBasecalled()) {
      return null;
    }

      // test if the fast5 file is barcoded
      if (isBarcoded()) {
          // test if the basecaller is Metrichor
          if (this.basecaller == Basecaller.METRICHOR) {

                  return reader.getStringAttribute("/Analyses/Barcoding_000/Barcoding",
                          "barcode_arrangement");

          }

          // test if the basecaller is Albacore
          if (this.basecaller == Basecaller.ALBACORE) {

              return reader.getStringAttribute("/Analyses/Barcoding_000/Summary/barcoding",
                      "barcode_full_arrangement");
          }
      }


    return null;
  }

  //
  //
  // FASTQ getters
  //
  //

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * template sequence in the fast5 file.
   * @return a string with the sequence fastq of the template strand
   */
  public String getTemplateFastq() {

    // test if the fast5 file is basecalled
    if (!isBasecalled()) {
      return null;
    }

    return fix(reader
            .readString("/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"));
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * complemente sequence in the fast5 file.
   * @return a string with the sequence fastq of the complemente strand
   */
  public String getComplementFastq() {

    // test if the fast5 file is basecalled and is 2D
    if (!is2D() || !isBasecalled()) {
      return null;
    }

    return fix(reader
            .readString("/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq"));
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * transcript sequence contains adaptor5'+transcript-consensus+adaptor3'.
   * @return a string with the sequence fastq of the transcript+rt-adaptor
   */
  public String getTranscriptFastq() {

    // test if the fast5 file is basecalled and is barcoded
    if (!isBarcoded() || !isBasecalled()) {
      return null;
    }

    return fix(reader.readString("/Analyses/Barcoding_000/Barcoding/Fastq"));
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * consensus sequence contains
   * barcodePos0+adaptor5'+transcript-consensus+adaptor3'+barcodePos1.
   * @return a string with the sequence fastq of the consensus
   */
  public String getConsensusFastq() {

    // test if the fast5 file is basecalled and is 2D
    if (!is2D() || !isBasecalled()) {
      return null;
    }

    return fix(
            reader.readString("/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"));
  }

  //
  //
  // Log Status getters
  //
  //

  /**
   * Method of the class Fast5 to obtain the status of the barcoding workflow.
   * @return a string of the status of the barcode workflow
   */
  public String getBarcodindFinalStatus() {

      // test if the fast5 file is basecalled and is barcoded
      if (!isBasecalled() || !isBarcoded()) {
          return null;
      }

      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR){
          return getLogFinalStatus(reader.readString("/Analyses/Barcoding_000/Log"));
      }

      //test if the basecaller is Albacore
      if(this.basecaller==Basecaller.ALBACORE){
          return getLogFinalStatus(reader.readString("/Analyses/Barcoding_000/log"));
      }

      return null;
  }

  /**
   * Method of the class Fast5 to obtain the status of the basecall1D workflow.
   * @return a string of the status of the basecall1D workflow
   */
  public String getBaseCall1DFinalStatus() {

    // test if the fast5 file is basecalled
    if (!isBasecalled()) {
      return null;
    }

      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          return getLogFinalStatus(
                  reader.readString("/Analyses/Basecall_1D_000/Log"));
      }
      return null;
  }

  /**
   * Method of the class Fast5 to obtain the status of the basecall2D workflow.
   * @return a string of the status of the basecall2D workflow
   */
  public String getBaseCall2DFinalStatus() {

    // test if the fast5 file is basecalled and is 2D
    if (!isBasecalled() || !is2D()) {
      return null;
    }
      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          return getLogFinalStatus(
                  reader.readString("/Analyses/Basecall_2D_000/Log"));
      }
      return null;
  }

  /**
   * Method of the class Fast5 to obtain the status of the Calibration Strand
   * workflow.
   * @return a string of the status of the Calibration Strand workflow
   */
  public String getCalibrationStrandFinalStatus() {

    // test if the fast5 file is basecalled
    if (!isBasecalled()) {
      return null;
    }

      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          return getLogFinalStatus(
                  reader.readString("/Analyses/Calibration_Strand_000/Log"));
      }

      return null;
  }

  /**
   * Method of the class Fast5 to obtain the status of the Event Detection
   * workflow.
   * @return a string of the status of the Event Detection workflow
   */
  public String getEventDetectionFinalStatus() {

    // test if the fast5 file is basecalled and the chemi is R7.3
    if (!isBasecalled() || getChemistryVersion() == ChemistryVersion.R7_3) {
      return null;
    }
      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          return getLogFinalStatus(
                  reader.readString("/Analyses/EventDetection_000/Log"));
      }
      return null;
  }

  /**
   * Method of the class Fast5 to obtain the status of the Hairpin split
   * workflow.
   * @return a string of the status of the Hairpin split workflow
   */
  public String getHairpinSplitFinalStatus() {

    // test if the fast5 file is basecalled and is 2D
    if (!isBasecalled() || !is2D()) {
      return null;
    }

      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          return getLogFinalStatus(
                  reader.readString("/Analyses/Hairpin_Split_000/Log"));
      }
      return null;
  }

  /**
   * Method of the class Fast5 to obtain the final status of the workflow.
   * @param log of a workflow
   * @return a string of the status of the Hairpin split workflow
   */
  private String getLogFinalStatus(String log) {

      // test if the basecaller is Metrichor
      if(this.basecaller==Basecaller.METRICHOR) {
          // split the log by '\n'
          String[] work = log.split("[\n]");

          // split the log to have the before last line
          String[] work2 = work[work.length - 2].split("\\s");
          String Status = "";

          // get the essential message of the log : the status
          for (int i = 2; i < work2.length; i++) {
              Status += work2[i] + " ";
          }

          // delete the end point of the message
          Status = Status.substring(0, Status.length() - 1);

          return Status;
      }
      return null;
  }

  /**
   * Method of the class Fast5 to fix the line break of fastq.
   * @param sequence a string sequence
   * @return a string with a "\n" at the end
   */
  private static String fix(String sequence) {

    // test if the sequence is null
    if (sequence == null) {
      return null;
    }

    // test if the sequence is equal to 1 in length
    if (sequence.length() == 1) {
      return "";
    }

    // return the sequence fastq corrected
    return sequence
            + (sequence.charAt(sequence.length() - 1) != '\n' ? "\n" : "");

  }

}