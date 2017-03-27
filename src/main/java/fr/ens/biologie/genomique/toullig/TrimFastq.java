package fr.ens.biologie.genomique.toullig;

import com.google.common.primitives.Ints;
import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import fr.ens.biologie.genomique.toullig.trimming.TrimModes;
import fr.ens.biologie.genomique.toullig.trimming.TrimWithCutadapt;
import fr.ens.biologie.genomique.toullig.trimming.TrimWithTrimmomatic;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;


import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static fr.ens.biologie.genomique.eoulsan.bio.Alphabets.*;
import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;

/**
 * Created by birer on 24/02/17.
 */
public class TrimFastq {

    private HashMap<String, String[]> fastqHash = new HashMap<String, String[]>();
    private String adaptorRT;
    private String adaptorStrandSwitching;
    private String pathOutputTrimLeftFasta;
    private String pathOutputTrimRightFasta;

    private File samFile;
    private File fastqFile;
    private File workDir;
    private File nameOutputFastq;
    private File adaptorFile;

    private boolean processPTrim=true;
    private boolean processSWTrim=false;
    private boolean processCutadapt=true;
    private boolean processTrimmomatic=false;
    private boolean processStats=false;

    private int addIndexOutlier=15;
    private int lengthWindowsSW=15;
    private double thresholdSW=0.8;
    private double errorRateCutadapt=0.4;
    private int seedMismatchesTrimmomatic=17;
    private int palindromeClipThresholdTrimmomatic=30;
    private int simpleClipThreshold=7;

    private BufferedWriter fastaFileLeftOutlier;
    private BufferedWriter fastaFileRightOutlier;


    public TrimFastq(File samFile, File fastqFile, File adaptorFile, File nameOutputFastq, File workDir) throws IOException, InterruptedException {

        if (!samFile.exists()) {
            throw new IOException(
                    "The file " + samFile + " dont exist!");
        } else {
            this.samFile = samFile;
        }

        if (!fastqFile.exists()) {
            throw new IOException(
                    "The file " + fastqFile + " dont exist!");
        } else {
            this.fastqFile = fastqFile;
        }

        if (!nameOutputFastq.exists()) {
            try{
                PrintWriter nameOutputFastqWriter = new PrintWriter(nameOutputFastq.toString(), "UTF-8");
                nameOutputFastqWriter.close();
            } catch (IOException e){
                e.printStackTrace();
            }

            this.nameOutputFastq = nameOutputFastq;
        } else {
            this.nameOutputFastq = nameOutputFastq;
        }

        if (!adaptorFile.exists()) {
            throw new IOException(
                    "The file " + adaptorFile + " dont exist!");
        } else {
            this.adaptorFile = adaptorFile;
        }

        if (!workDir.exists()) {
            throw new IOException(
                    "The directory " + workDir + " dont exist!");
        } else {
            this.workDir = workDir;
        }
    }

    //
    // Read Files
    //

    /**
     * Method of the class TrimFastq to read and get the adaptor in the adaptor file.
     * @param adaptorFile, the file who contains adaptor
     * @throws IOException
     */
    private void readAdaptorRTFile(BufferedReader adaptorFile)throws IOException{

        int i =0;
        String line;
        while ((line = adaptorFile.readLine()) != null) {
            if(i==1){
                this.adaptorRT =line;
            }
            if(i==3){
                this.adaptorStrandSwitching =line;
            }
            i++;
        }
        adaptorFile.close();
    }

    /**
     * Method of the class TrimFastq to read and get the ID,sequence, score and CIGAR of a sam file.
     * @param samInputStream, the input stream file sam
     * @throws IOException
     */
    private void readSamFile(InputStream samInputStream) throws IOException {

        try(final SamReader inputSam =
                    SamReaderFactory.makeDefault().open(SamInputResource.of(samInputStream))){

            String lengthBeginOutlier="";
            String lengthEndOutlier="";
            String score="";
            String sequence="";

            HashMap<String, String> fastqHashMultimapped = new HashMap<String, String>();

            for (SAMRecord samRecord : inputSam) {

                String CIGAR=samRecord.getCigarString();

                String QFlag = ""+samRecord.getFlags();
                String ID =samRecord.getReadName();
                int cigarLength=samRecord.getCigarLength();


                if(this.fastqHash.containsKey(ID)) {
                    fastqHashMultimapped.put(ID,"");
                    String[] tabvalue= this.fastqHash.get(ID);
                    if(Integer.parseInt(tabvalue[6])<=cigarLength){
                        this.fastqHash.put(ID, new String[]{sequence, score, CIGAR, lengthBeginOutlier, lengthEndOutlier, QFlag, ""+cigarLength});
                    }else{
                        continue;
                    }

                }else {
                    this.fastqHash.put(ID, new String[]{sequence, score, CIGAR, lengthBeginOutlier, lengthEndOutlier, QFlag, ""+cigarLength});
                }
            }
            System.out.println("Nombre de sequence qui ont fait ont mappé plusieurs fois : "+fastqHashMultimapped.size());
            inputSam.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }


    /**
     * Method of the class TrimFastq to read a fastq file.
     * @param fastqBufferedReader, a BufferedReader fastq file
     * @throws IOException
     */
    private void readFastqFile(File fastqBufferedReader)throws IOException{

        try (FastqReader reader = new FastqReader(fastqBufferedReader)) {
            for (ReadSequence read : reader) {

                String header = read.getName();
                String[] part=header.split(" ");
                String ID = part[0];
                String sequence = read.getSequence();
                String score = read.getQuality();
                String[] tabValue = this.fastqHash.get(ID);
                tabValue[0]=sequence;
                tabValue[1]=score;
            }
            reader.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    //
    //  Set
    //

    /**
     * Method of the class TrimFastq to set the process P in trim mode.
     * @param processPTrim, a boolean
     */
    public void setProcessPTrim(boolean processPTrim) {
        this.processPTrim = processPTrim;
    }

    /**
     * Method of the class TrimFastq to set the process SW in trim mode.
     * @param processSWTrim, a boolean
     */
    public void setProcessSWTrim(boolean processSWTrim) {
        this.processSWTrim = processSWTrim;
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
     * @param processTrimmomatic, a boolean
     */
    public void setProcessTrimmomatic(boolean processTrimmomatic) {
        this.processTrimmomatic = processTrimmomatic;
    }

    /**
     * Method of the class TrimFastq to set the process stats.
     * @param processStats, a boolean
     */
    public void setProcessStats(boolean processStats) {
        this.processStats = processStats;
    }

    public void setLengthWindowsSW(int lengthWindowsSW){
        this.lengthWindowsSW=lengthWindowsSW;
    }

    public void setThresholdSW(double thresholdSW){
        this.thresholdSW=thresholdSW;
    }

    public void setErrorRateCutadapt(double errorRateCutadapt){
        this.errorRateCutadapt=errorRateCutadapt;
    }

    public void setSeedMismatchesTrimmomatic(int seedMismatchesTrimmomatic){
        this.seedMismatchesTrimmomatic=seedMismatchesTrimmomatic;
    }

    public void setPalindromeClipThresholdTrimmomatic(int palindromeClipThresholdTrimmomatic){
        this.palindromeClipThresholdTrimmomatic=palindromeClipThresholdTrimmomatic;
    }

    public void setSimpleClipThreshold(int simpleClipThreshold){
        this.simpleClipThreshold=simpleClipThreshold;
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
    private void executionTrimWithCutadapt(File fastqFile) throws IOException, InterruptedException {

        String pathFastaFileLeftOutlier = this.workDir+"/fastaFileLeftOutlier.fasta";
        String pathFastaFileRightOutlier = this.workDir+"/fastaFileRightOutlier.fasta";

        String infoTrimLeftPath = this.workDir+"/logCutadaptLeftOutlier.txt";
        String infoTrimRightPath = this.workDir+"/logCutadaptRightOutlier.txt";

        readFastqFile(fastqFile);

        this.fastaFileLeftOutlier= new BufferedWriter(new FileWriter(pathFastaFileLeftOutlier));
        this.fastaFileRightOutlier= new BufferedWriter(new FileWriter(pathFastaFileRightOutlier));

        System.out.println("Trim begin !");

        TrimModes trimMode = new TrimModes(this.lengthWindowsSW, this.thresholdSW, this.fastqHash, this.nameOutputFastq, this.addIndexOutlier,this.adaptorFile, this.seedMismatchesTrimmomatic, this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold, this.processCutadapt, this.processTrimmomatic, pathFastaFileLeftOutlier, pathFastaFileRightOutlier);

        if(this.processPTrim){
            this.fastqHash=trimMode.trimOutlierWithPMethod();
        }
        if(this.processSWTrim){
            this.fastqHash=trimMode.trimOutlierWithSWMethod();
        }

        this.fastaFileRightOutlier.close();
        this.fastaFileLeftOutlier.close();

        TrimWithCutadapt trimmingCutadapt = new TrimWithCutadapt(this.fastqHash, this.nameOutputFastq, this.pathOutputTrimLeftFasta, this.pathOutputTrimRightFasta, this.adaptorRT, this.adaptorStrandSwitching, this.errorRateCutadapt);

        System.out.println("Begin use cutadapt !");

        String strandLeft = "-g";
        String strandRight = "-a";

        // Cutadapt execution
        trimmingCutadapt.cutadaptTrim(pathFastaFileLeftOutlier, strandLeft, infoTrimLeftPath, this.pathOutputTrimLeftFasta);
        trimmingCutadapt.cutadaptTrim(pathFastaFileRightOutlier, strandRight, infoTrimRightPath, this.pathOutputTrimRightFasta);

        // Merge the output form cutadapt
        trimmingCutadapt.mergeTrimOutlier();

        if (processStats) {
            System.out.println("Start stat Left outlier");
            trimmingCutadapt.statsLogCutadapt(infoTrimLeftPath);
            System.out.println("Start stat Right outlier");
            trimmingCutadapt.statsLogCutadapt(infoTrimRightPath);
        }

    }

    /**
     * Method of the class TrimFastq to execute the process with trimmomatic.
     * @param fastqFile, a fastq file
     * @throws IOException
     * @throws InterruptedException
     */
    private void executionTrimWithTrimmomatic(File fastqFile) throws IOException, InterruptedException {

        readFastqFile(fastqFile);

        System.out.println("Begin use trimmomatic !");

        TrimModes trimMode = new TrimModes(this.lengthWindowsSW, this.thresholdSW, this.fastqHash, this.nameOutputFastq, this.addIndexOutlier,this.adaptorFile, this.seedMismatchesTrimmomatic, this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold, this.processCutadapt, this.processTrimmomatic, "","");

        if(this.processPTrim){
            this.fastqHash=trimMode.trimOutlierWithPMethod();
        }
        if(this.processSWTrim){
            this.fastqHash=trimMode.trimOutlierWithSWMethod();
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

        this.pathOutputTrimLeftFasta= this.workDir+"/outputFastaFileLeftOutlier.fastq";
        this.pathOutputTrimRightFasta= this.workDir+"/outputFastaFileRightOutlier.fastq";

        //Problem with ONT skip read for the RT adaptor (to many TTTT..)
        readAdaptorRTFile(new BufferedReader(new FileReader(this.adaptorFile)));

        InputStream samInputStream =new FileInputStream(this.samFile);
        readSamFile(samInputStream);

        if(this.processCutadapt){
            executionTrimWithCutadapt(this.fastqFile);
        }
        if(this.processTrimmomatic){
            executionTrimWithTrimmomatic(this.fastqFile);
        }
    }
}

