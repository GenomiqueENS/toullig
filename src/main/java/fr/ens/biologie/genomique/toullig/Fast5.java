package fr.ens.biologie.genomique.toullig;

import java.io.File;
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
  private final IHDF5Reader reader;

  /**
   * Constructor of the Fast5 class.
   * @param fast5File constructor
   */
  Fast5(File fast5File) {

    this.reader = readFast5File(fast5File);
    this.status = readStatus();
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
  };

  /**
   * Values of the variable Type that design the type of experimental design.
   */
  public enum Type {
    TYPE_1D, TYPE_2D // Type of the sequencing
  };

  /**
   * Values of the variable Status that design the state of the fast5 file (if
   * it was basecalled or not).
   */
  public enum Status {
    PRE_BASECALLING, AFTER_BASECALLING // The status of the fast5 file
  };

  /**
   * Values of the variable RVersion that design the chemical kit use.
   */
  public enum ChemistryVersion {
    R7_3, R9, R9_4 // Chemical version of the kit
  };

  //
  // read files
  //

  /**
   * Method who use the class IDFH5Reader to read the fast5 file.
   * @return a hdf5 file open
   */

  private static IHDF5Reader readFast5File(File fast5File) {
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
    if (reader
            .isGroup("/Analyses/Basecall_2D_000/BaseCalled_template") == true) {
      return Version.V1_0;
    }
    if (reader
            .isGroup("/Analyses/Basecall_1D_000/BaseCalled_template") == true) {
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
    if (reader.isGroup("/Analyses/Basecall_2D_000") == true) {
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
    if (this.reader.getFile() == null) {
      throw new IllegalStateException("The file is closed");
    }
    if (reader.isGroup("/Analyses/Basecall_1D_000") == true) {
      return Status.AFTER_BASECALLING;
    } else {
      return Status.PRE_BASECALLING;
    }
  }

  /**
   * Method of the variable RVersion who obtain the chemical kit use.
   * @return a rversion with the chemical kit version use
   */

  private ChemistryVersion readChemistryVersion() {

    // Case of not basecalled fast5 file
    if (!isBasecalled() && reader.isGroup("/Raw/Reads")) {
      return ChemistryVersion.R9;
    }
    if (!isBasecalled()
            && reader.isGroup("/Analyses/EventDetection_000/Reads")) {
      return ChemistryVersion.R7_3;
    }

    // Case of basecalled fast5 file

    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r7.3_")) {
      return ChemistryVersion.R7_3;
    }

    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r9_")) {
      return ChemistryVersion.R9;

    }
    if (reader
            .getStringAttribute("/Analyses/Basecall_1D_000/Configuration/general",
                    "model_type")
            .contains("r94_")) {
      return ChemistryVersion.R9_4;

    }
    return null;
  }

  //
  //
  // Important getters
  //
  //

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
    return reader.isGroup("/Analyses/Barcoding_000");// Si l'échantillion a
    // barcodé
  }

  /**
   * Boolean shortcut to know if the file is basecalled.
   * @return a boolean with the basecalled information
   */
  public boolean isBasecalled() {
    return this.status == Status.AFTER_BASECALLING;// Si l'échantillion a
    // été analysé
  }

  /**
   * Boolean shortcut to know if the experiement is 2D (or 1D in the opposite
   * case).
   * @return a boolean with the type 2D information
   */
  public boolean is2D() {
    return this.type == Type.TYPE_2D;// Si l'échantillion a
    // été analysé
  }

  //
  //
  // primary information getters
  //
  //

  //
  // Raw Group Information
  //

  /**
   * Method of the class Fast5 to obtain the electrical signal data of the fast5
   * file.
   * @return a double tabular with the electrical signal
   */
  public int[] getElectricalSignal() {
    if (this.chemistryVersion == ChemistryVersion.R9) {
      return reader
              .readIntArray("/Raw/Reads/Read_" + getNumberRead() + "/Signal");
    }
    return null;
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
  public Date getdDateExp() {
    String dateInt = reader.getStringAttribute("/UniqueGlobalKey/tracking_id",
            "exp_start_time");
    return new Date(Long.parseLong(dateInt) * 1000);
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
    if (!isBasecalled() && getChemistryVersion() == ChemistryVersion.R9) {
      String s = reader.getAllGroupMembers("/Raw/Reads").get(0);
      return Integer.parseInt(s.substring(s.indexOf('_') + 1));
    }
    if (!isBasecalled() && getChemistryVersion() == ChemistryVersion.R7_3) {
      String s = reader.getAllGroupMembers("/Analyses/EventDetection_000/Reads")
              .get(0);
      return Integer.parseInt(s.substring(s.indexOf('_') + 1));
    }
    return Integer.parseInt(reader.getStringAttribute(
            "/Analyses/Basecall_1D_000/Configuration/general", "read_id"));
  }

  /**
   * Method of the class Fast5 to obtain the name of the basecaller and the
   * version in the fast5 file.
   * @return a string with the basecaller name and version
   */
  public String getBaseCaller() {
    if (!isBasecalled()) {
      return null;
    }
    String chimaeraVersion = reader
            .getStringAttribute("/Analyses/Basecall_1D_000", "chimaera version");
    String dragonetVersion = reader
            .getStringAttribute("/Analyses/Basecall_1D_000", "dragonet version");
    return "chimaera v" + chimaeraVersion + " | dragonet v" + dragonetVersion;
  }

  /**
   * Method of the class Fast5 to obtain alignement of the complement and the
   * template sequence in the fast5 file.
   * @return a string with the alignement
   */
  public String getAlignment() {
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
    if (!isBasecalled()) {
      return null;
    }
    if (isBarcoded()) {
      return reader.getStringAttribute("/Analyses/Barcoding_000/Barcoding",
              "barcode_arrangement");
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
    if (!is2D() || !isBasecalled()) {
      return null;
    }
    return fix(reader.readString("/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq"));
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * transcript sequence contains adaptor5'+transcript-consensus+adaptor3'.
   * @return a string with the sequence fastq of the transcript+rt-adaptor
   */
  public String getTranscriptFastq() {
    if (!isBarcoded() || !isBasecalled()) {
      return null;
    }
    return fix(reader.readString("/Analyses/Barcoding_000/Barcoding/Fastq"));
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * consensus sequence contains barcodePos0+adaptor5'+transcript-consensus+adaptor3'+barcodePos1.
   * @return a string with the sequence fastq of the consensus
   */
  public String getConsensusFastq() {
    if (!is2D() || !isBasecalled()) {
      return null;
    }
    return fix(reader.readString("/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"));
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
    if (!isBasecalled() || !isBarcoded()) {
      return null;
    }
    return getLogFinalStatus(reader.readString("/Analyses/Barcoding_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the status of the basecall1D workflow.
   * @return a string of the status of the basecall1D workflow
   */
  public String getBaseCall1DFinalStatus() {
    if (!isBasecalled()) {
      return null;
    }
    return getLogFinalStatus(
            reader.readString("/Analyses/Basecall_1D_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the status of the basecall2D workflow.
   * @return a string of the status of the basecall2D workflow
   */
  public String getBaseCall2DFinalStatus() {
    if (!isBasecalled() || !is2D()) {
      return null;
    }
    return getLogFinalStatus(
            reader.readString("/Analyses/Basecall_2D_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the status of the Calibration Strand
   * workflow.
   * @return a string of the status of the Calibration Strand workflow
   */
  public String getCalibrationStrandFinalStatus() {
    if (!isBasecalled()) {
      return null;
    }
    return getLogFinalStatus(
            reader.readString("/Analyses/Calibration_Strand_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the status of the Event Detection
   * workflow.
   * @return a string of the status of the Event Detection workflow
   */
  public String getEventDetectionFinalStatus() {
    if (!isBasecalled() || getChemistryVersion() == ChemistryVersion.R7_3) {
      return null;
    }
    return getLogFinalStatus(
            reader.readString("/Analyses/EventDetection_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the status of the Hairpin split
   * workflow.
   * @return a string of the status of the Hairpin split workflow
   */
  public String getHairpinSplitFinalStatus() {
    if (!isBasecalled() || !is2D()) {
      return null;
    }
    return getLogFinalStatus(
            reader.readString("/Analyses/Hairpin_Split_000/Log"));
  }

  /**
   * Method of the class Fast5 to obtain the final status of the workflow.
   * @param log of a workflow
   * @return a string of the status of the Hairpin split workflow
   */
  public String getLogFinalStatus(String log) {
    String[] work = log.split("[\n]");
    String[] work2 = work[work.length - 2].split("\\s");
    String Status = "";
    for (int i = 2; i < work2.length; i++) {
      Status += work2[i] + " ";
    }
    Status = Status.substring(0, Status.length() - 1);
    return Status;
  }

  /**
   * Method of the class Fast5 to fix the line break of fastq.
   * @param s a string
   * @return a string with a "\n" at the end
   */
  private static String fix(String s) {

    if (s==null) {
      return null;
    }

    if (s.length()==1) {
      return "";
    }

    return s + (s.charAt(s.length()-1) != '\n' ? "\n" : "");

  }

}