package fr.ens.biologie.genomique.toullig.GtftoGdp;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import fr.ens.biologie.genomique.eoulsan.bio.GFFEntry;
import fr.ens.biologie.genomique.eoulsan.bio.io.GFFReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.GTFReader;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

/**
 * Created by birer on 20/06/17. Class to translate GTF file format to GPD file
 * format
 */
public class GtftoGpd {

  private File gtfFile;
  private File gpdfile;

  public GtftoGpd(File gtfFile, File gpdfile) throws IOException {

      // test if the GTF file exist
      if(gtfFile.exists()){
          this.gtfFile = gtfFile;
      }

      this.gpdfile=gpdfile;

  }

    /**
     * Method of the class GtftoGpd who obtain the data to write a GPD format in a GTF file format.
     * @return HashMap<String, List> contains all data for the GPD format
     * @throws FileNotFoundException to test if the gtf file exist.
     */
  public HashMap<String, List> getGtfFile() throws FileNotFoundException {

    InputStream annotationIs = new FileInputStream(this.gtfFile);

    boolean gtfFormat = true;

    HashMap<String, List> workHash = new HashMap<>();

    final GFFReader gffReader =
        gtfFormat ? new GTFReader(annotationIs) : new GFFReader(annotationIs);

    String pivotTranscriptId = "";
    int firstLoop = 0;

    int countExon = 0;

    List<Integer> startPositionExonList = new ArrayList<>();
    List<Integer> endPositionExonList = new ArrayList<>();

    int counter = 0;

    String type = "";
    String chrom = "";
    int start = 0;
    int end = 0;
    String strand = "";
    String geneId = "";
    String transcriptID = "";

    for (final GFFEntry gff : gffReader) {

      // System.out.println(gff);
      counter++;

      String part[] = gff.toString().split("\t");
      type = part[2];

      String part2[] = part[8].split(";");
      if (part2[2].contains("transcript_id")) {
        transcriptID = part2[2].split("=")[1];
      }

      if (type.equals("transcript")) {

        // Initialize for the first loop
        if (firstLoop == 0) {

          pivotTranscriptId = transcriptID;
          firstLoop = 1;
        }
      }

      // treat a transcript
      if (pivotTranscriptId.equals(transcriptID)) {

        if (type.equals("transcript")) {

          chrom = part[0];
          start = Integer.parseInt(part[3]);
          end = Integer.parseInt(part[4]);
          strand = part[6];
          geneId = part2[0].split("=")[1];
        }

        if (type.equals("exon")) {

          start = Integer.parseInt(part[3]);
          end = Integer.parseInt(part[4]);

          countExon += 1;
          startPositionExonList.add(start);
          endPositionExonList.add(end);

        }
      }
      // new transcript re-initialize + workhash gathering
      else {

        List<String> informationTranscript = new ArrayList<String>();
        informationTranscript.add(geneId
            + ";" + chrom + ";" + strand + ";" + start + ";" + end + ";"
            + countExon + ";" + startPositionExonList + ";"
            + endPositionExonList);
        workHash.put(pivotTranscriptId, informationTranscript);

        // reset variables
        countExon = 0;
        startPositionExonList.clear();
        endPositionExonList.clear();

        // initialize for the new transcript
        if (type.equals("transcript")) {

          chrom = part[0];
          start = Integer.parseInt(part[3]);
          end = Integer.parseInt(part[4]);
          strand = part[6];
          geneId = part2[0].split("=")[1];
        }

      }

      // re-initialize the 'pivot'
      pivotTranscriptId = transcriptID;

    }

    // For the last transcript
    List<String> informationTranscript = new ArrayList<String>();
    informationTranscript.add(geneId
        + ";" + chrom + ";" + strand + ";" + start + ";" + end + ";" + countExon
        + ";" + startPositionExonList + ";" + endPositionExonList);
    workHash.put(pivotTranscriptId, informationTranscript);

    // return hash
    return workHash;
  }

  //
  // Write
  //

    /**
     * Method of the class GtftoGpd who write the GTF format file.
     * @param workHash, the hash who contains all GPD data
     * @param gpdOutputFile, the output File
     * @throws IOException, if an error occur
     */
  public void writeGPD(HashMap workHash)
      throws IOException {

    Writer gpdWriter = new FileWriter(this.gpdfile);

    // read the work hash
    for (Object transcriptId : workHash.keySet()) {

      // get the informations transcript
      List<String> informationTranscript =
          (List<String>) workHash.get(transcriptId);
      String information = informationTranscript.toString().replace("[", "");
      information = information.replace("]", "");
      information = information.replace(", ", ",");
      String part[] = information.split(";");
      String geneId = part[0];
      String chrom = part[1];
      String strand = part[2];
      String start = part[3];
      String end = part[4];
      String countExon = part[5];
      String startPositionExonList = part[6];
      String endPositionExonList = part[7];

      // write the transcript line in gpd format
      gpdWriter.write(transcriptId
          + "\t" + geneId + "\t" + chrom + "\t" + strand + "\t" + start + "\t"
          + end + "\t" + start + "\t" + end + "\t" + countExon + "\t"
          + startPositionExonList + "\t" + endPositionExonList + "\n");
    }

    gpdWriter.close();

  }

  //
  // Main
  //

  public void execution() throws IOException {

    getLogger().info("GTF to GPD translation in progress...");

    // get the Workhash with all the data for the gpd format
    HashMap<String, List> workHash = getGtfFile();

    // write the gpd format file
    writeGPD(workHash);

  }
}
