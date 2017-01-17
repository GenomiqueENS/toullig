package fr.ens.biologie.genomique.nanopore;

import java.io.File;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 * @author Aurelien Birer
 */
public class Fast5 implements AutoCloseable {

  private static final Pattern PATTERN1 = Pattern.compile("chimaera .*");
  private static final Pattern PATTERN2 = Pattern.compile("dragonet .*");
  private final Version version;
  private Type type;
  private final Status status;
  private RVersion rversion;
  private IHDF5Reader reader;

  /**
   * Constructor of the Fast5 class.
   * @param fast5File constructor
   */
  Fast5(File fast5File) {

    this.reader = readFast5File(fast5File);
    this.status = readStatus();
    this.version = readVersion();
    this.type = readType();
    this.rversion = readRVersion();

  }

  //
  // Principals declaration
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
  public enum RVersion {
    R7, R9, R9_4 // Chemical version of the kit
  };

  //
  // read files
  //

  /**
   * Method who use the class IDFH5Reader to read the fast5 file.
   * @return a hdf5 file open
   */

  public IHDF5Reader readFast5File(File fast5File) {
    IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
    return this.reader = hdf5Fac.openForReading(fast5File);
  }

  /**
   * Method of the variable Version who obtain the version of the fast5 format.
   * @return a version with the version of the fast5 format
   */
  public Version readVersion() {
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
  public Type readType() {
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
  public Status readStatus() {
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
  public RVersion readRVersion() {
    if (!isBasecalled() && reader.isGroup("/Raws/Reads")) {
      return RVersion.R9;
    }
    if (!isBasecalled()
        && reader.isGroup("/Analyses/EventDetection_000/Reads")) {
      return RVersion.R7;
    }
    if (reader.isGroup("/Analyses/Alignment_000") == true) {
      return RVersion.R7;
    }
    if (reader.isGroup(
        "/Analyses/Basecall_1D_000/Configuration/genome_mapping") == false
        && reader.isGroup("/Analyses/Alignment_000") == false) {
      return RVersion.R9;
    }
    if (reader.isGroup(
        "/Analyses/Basecall_1D_000/Configuration/genome_mapping") == true) {
      return RVersion.R9_4;
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
  public RVersion getRversion() {
    return this.rversion;
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
    // if (this.rversion == RVersion.R7) {
    // System.out.println(
    // reader.isDataSet("/Analyses/EventDetection_000/Reads/Read_" +
    // getNumberRead() + "/Events"));
    // //
    // System.out.println(reader.("/Analyses/EventDetection_000/Reads/Read_"
    // // + getNumberRead() + "/Events"));
    // return null;
    // }
    if (this.rversion == RVersion.R9) {
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
   * Method of the class Fast5 to obtain the flowcell version in the fast5 file.
   * @return a string with the flowcell version
   */
  public String getFlowcellVersion() {
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
    if (this.rversion == RVersion.R7 || !isBasecalled()) {
      return null;

    }
    return reader.getStringAttribute("/UniqueGlobalKey/context_tags",
        "experiment_type");
  }

  /**
   * Method of the class Fast5 to obtain the frequency of the sample in the
   * fast5 file.
   * @return an int with sample frequency
   */
  public int getSampleFrequency() {
    return Integer.parseInt(reader.getStringAttribute("/UniqueGlobalKey/context_tags",
        "sample_frequency"));
  }

  /**
   * Method of the class Fast5 to obtain the channel number of the pore in the
   * fast5 file.
   * @return an int with sample frequency
   */
  public int getChannelNumber() {
    return Integer.parseInt(reader.getStringAttribute("/UniqueGlobalKey/channel_id",
        "channel_number"));
  }

  //
  //
  // basecalling information getters
  //
  //

  /**
   * Method of the class Fast5 to obtain the number of the read in the fast5
   * file.
   * @return an int with number of the read
   */
  public int getNumberRead() {
    if (!isBasecalled() && getRversion() == RVersion.R9) {
      String s = reader.getAllGroupMembers("/Raw/Reads").get(0);
      return Integer.parseInt(s.substring(s.indexOf('_')+1));
      
    }
    if (!isBasecalled() && getRversion() == RVersion.R7) {
      String s = reader.getAllGroupMembers("/Analyses/EventDetection_000/Reads").get(0);
      return Integer.parseInt(s.substring(s.indexOf('_')+1));
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
    String log = reader.readString("/Analyses/Basecall_1D_000/Log");
    String Basecallers = "";
    Matcher matcher1 = PATTERN1.matcher(log);
    while (matcher1.find()) {
      Basecallers = matcher1.group();
    }

    Matcher matcher2 = PATTERN2.matcher(log);
    while (matcher2.find()) {
      Basecallers = Basecallers + " | " + matcher2.group();
    }
    return Basecallers;
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
    return reader
        .readString("/Analyses/Basecall_1D_000/BaseCalled_template/Fastq");
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
    return reader
        .readString("/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq");
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * hairpin sequence in the fast5 file.
   * @return a string with the sequence fastq of the hairpin
   */
  public String getHairpinFastq() {
    if (!is2D() || !isBasecalled()) {
      return null;
    }
    // If the file is fail after the basecalling online, it's possible that the
    // group "/Analyses/Basecall_2D_000/BaseCalled_2D" dont exist.
    return reader.readString("/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq");
  }

  /**
   * Method of the class Fast5 to obtain the sequence fastq + score of the
   * barcode sequence in the fast5 file.
   * @return a string with the sequence fastq of the barcode
   */
  public String getBarcodingFastq() {
    if (!isBarcoded() || !isBasecalled()) {
      return null;
    }
    return reader.readString("/Analyses/Barcoding_000/Barcoding/Fastq");
  }

}