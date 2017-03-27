package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import fr.ens.biologie.genomique.toullig.Utils;

import java.io.*;
import java.util.HashMap;
import java.util.concurrent.TimeUnit;

import static fr.ens.biologie.genomique.eoulsan.bio.Alphabets.AMBIGUOUS_DNA_ALPHABET;
import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;

/**
 * Created by birer on 27/03/17.
 */
public class TrimWithCutadapt {

    private HashMap<String, String[]> fastqHash;
    private File nameOutputFastq;
    private String pathOutputTrimLeftFasta;
    private String pathOutputTrimRightFasta;
    private String adaptorRT;
    private String adaptorStrandSwitching;
    private Alphabet alphabet = AMBIGUOUS_DNA_ALPHABET;
    private double errorRateCutadapt;

    public TrimWithCutadapt(HashMap fastqHash, File nameOutputFastq, String pathOutputTrimLeftFasta, String pathOutputTrimRightFasta, String adaptorRT, String adaptorStrandSwitching, double errorRateCutadapt) {

        this.fastqHash=fastqHash;
        this.nameOutputFastq=nameOutputFastq;
        this.pathOutputTrimLeftFasta=pathOutputTrimLeftFasta;
        this.pathOutputTrimRightFasta=pathOutputTrimRightFasta;
        this.adaptorRT=adaptorRT;
        this.adaptorStrandSwitching=adaptorStrandSwitching;
        this.errorRateCutadapt=errorRateCutadapt;

    }

    //
    //Merge outlier
    //

    /**
     * Method of the class TrimFastq to merge the results of cutadapt with the trim sequence.
     * @throws IOException
     */
    public void mergeTrimOutlier() throws IOException {

        BufferedWriter fastqTrimBufferedWritter= new BufferedWriter(new FileWriter(this.nameOutputFastq));
        HashMap<String, String> fastaLeftHash = new HashMap<String, String>();
        HashMap<String, String> fastaRightHash = new HashMap<String, String>();
        String shortestFastqSequence="";
        int i = 0;

        //Read of the left outlier fasta output by cutadapt
        try{
            BufferedReader outputTrimLeftFastaReader = new BufferedReader(new FileReader(this.pathOutputTrimLeftFasta));
            String line = "";
            while ((line = outputTrimLeftFastaReader.readLine()) != null) {
                if(line.contains(">")){
                    String[] part=line.split(">");
                    fastaLeftHash.put(part[1],outputTrimLeftFastaReader.readLine());
                }
            }
        }catch (IOException e){
            e.printStackTrace();
        }

        //Read of the right outlier fasta output by cutadapt
        try{
            BufferedReader outputTrimRightFastaFile = new BufferedReader(new FileReader(this.pathOutputTrimRightFasta));
            String line = "";
            while ((line = outputTrimRightFastaFile.readLine()) != null) {
                if(line.contains(">")){
                    String[] part=line.split(">");
                    fastaRightHash.put(part[1],outputTrimRightFastaFile.readLine());
                }
            }
        }catch (IOException e){
            e.printStackTrace();
        }

        System.out.println("Start of merging sequence !");

        //merge of the outlier trim on the main sequence
        for (String ID : this.fastqHash.keySet()) {

            i++;
            String leftSequence="";
            String rightSequence="";
            int lengthBeginOutlier=0;
            int lengthEndOutlier=0;
            String[] tabValue = this.fastqHash.get(ID);
            String cigar=tabValue[2];
            String score=tabValue[1];
            String sequence = tabValue[0];
            String beginOutlier= tabValue[3];
            String endOutlier= tabValue[4];

            if(!cigar.equals("*")) {

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

                if(lengthBeginOutlier<0){
                    lengthBeginOutlier=0;
                }
                if(lengthEndOutlier<0){
                    lengthEndOutlier=0;
                }

                String mainSequenceWithoutOutlier ="";

                if((lengthBeginOutlier+lengthEndOutlier)<sequence.length()){
                    mainSequenceWithoutOutlier = sequence.substring(lengthBeginOutlier, sequence.length() - lengthEndOutlier);
                }

                if (fastaLeftHash.get(ID) == null) {
                    leftSequence = "";
                } else {
                    leftSequence = fastaLeftHash.get(ID);
                }
                if (fastaLeftHash.get(ID) == null) {
                    rightSequence = "";
                } else {
                    rightSequence = fastaRightHash.get(ID);
                }

                if(rightSequence == null){
                    rightSequence="";
                }

                if(leftSequence == null){
                    leftSequence="";
                }

                int leftLength=lengthBeginOutlier-leftSequence.length();
                int rightLength=(sequence.length()-lengthEndOutlier)+rightSequence.length();

                if(lengthBeginOutlier-leftSequence.length()<0){
                    leftLength=0;
                }

                if(rightLength>score.length()){
                    rightLength=score.length();
                }

                // case with the addition of index (>0)
                if(leftLength>rightLength){
                    leftLength=rightLength;
                }

                String sequenceTrim = sequence.substring(leftLength,rightLength);

                String scoreTrim = score.substring(leftLength,rightLength );

                if(scoreTrim.length()!=sequenceTrim.length()){
                    System.out.println("problem :  "+scoreTrim.length()+"     "+sequenceTrim.length()+"   "+mainSequenceWithoutOutlier.length());
                    System.out.println(leftLength+"     "+rightLength+"     "+ID);
                }

                if(!sequenceTrim.equals("")){
                    Utils util = new Utils();
                    util.writeFastq(ID, sequenceTrim, scoreTrim, fastqTrimBufferedWritter);
                }

                if (i == 1) {
                    shortestFastqSequence = sequenceTrim;
                }
                if (shortestFastqSequence.length() >= sequenceTrim.length()) {
                    shortestFastqSequence = sequenceTrim;
                }
            }
        }
        System.out.println("la séquence trim la plus courte fait : " + shortestFastqSequence.length());
        fastqTrimBufferedWritter.close();
    }

    //
    //Cutadapt
    //

    /**
     * Method of the class TrimFastq to execute cutadapt.
     * @param pathFastaFileOutlier, the path of the fasta file with outlier to trim
     * @param strand, the cutadapt strand of the adaptor
     * @param infoTrimPath, path to the output log of cutadapt
     * @throws IOException
     * @throws InterruptedException
     */
    public void cutadaptTrim(String pathFastaFileOutlier, String strand, String infoTrimPath, String pathOutputTrimLeftFasta) throws IOException, InterruptedException {

        try {
            Utils utils = new Utils();
            String reverseComplementAdaptorRT=strand+" reverse_complement_RT_adaptor="+reverseComplement(this.adaptorRT,this.alphabet);
            String reverseComplementAdaptorSwithStrand=strand+" reverse_complement_Switch_Strand_RT_adaptor="+reverseComplement(this.adaptorStrandSwitching,this.alphabet);

            String complementAdaptorRT=strand+" complement_RT_adaptor="+ Utils.complement(this.adaptorRT,this.alphabet);
            String complementAdaptorSwithStrand=strand+" complement_Switch_Strand_RT_adaptor="+utils.complement(this.adaptorStrandSwitching,this.alphabet);

            StringBuffer reverseAdaptorRTStringBuffer = new StringBuffer(this.adaptorRT);
            StringBuffer reverseAdaptorSwithStrandStringBuffer = new StringBuffer(this.adaptorStrandSwitching);
            String reverseAdaptorRT= strand+" reverse_RT_adaptor="+reverseAdaptorRTStringBuffer.reverse().toString();
            String reverseAdaptorSwithStrand= strand+" reverse_Switch_Strand_RT_adaptor="+reverseAdaptorSwithStrandStringBuffer.reverse().toString();

            String adaptorRT=strand+" RT_adaptor="+this.adaptorRT;
            String adaptorStrandSwitching=strand+" Switch_Strand_RT_adaptor="+this.adaptorStrandSwitching;

            String quiet="";
            //String quiet="--quiet";

            ProcessBuilder pb = new ProcessBuilder("/bin/bash", "-c", "cutadapt "+adaptorRT+" "+adaptorStrandSwitching+" "+reverseAdaptorRT+" "+reverseAdaptorSwithStrand+" "+complementAdaptorRT+" "+complementAdaptorSwithStrand+" "+reverseComplementAdaptorRT+" "+reverseComplementAdaptorSwithStrand+" "+quiet+" --error-rate="+this.errorRateCutadapt+" --info-file="+infoTrimPath+" --overlap=7 --times=8 --match-read-wildcards --format=fasta "+pathFastaFileOutlier+" > "+pathOutputTrimLeftFasta);

            //System.out.println(pb.command());
            pb.redirectErrorStream(true);
            Process proc = pb.start(); // Start the process.

            getLogCutadapt(proc);

            // int exitproc=proc.waitFor();
            if(!proc.waitFor(1, TimeUnit.MINUTES)) {
                //timeout - kill the process.
                proc.destroy(); // consider using destroyForcibly instead
            }
        }
        catch (IOException e) {
            e.printStackTrace(); // or log it, or otherwise handle it
        }
        catch (InterruptedException ie) {
            ie.printStackTrace(); // or log it, or otherwise handle it
        }
    }

    /**
     * Method of the class TrimFastq to display log of cutadapt (delete the option --quiet).
     * @param proc, a Processus
     * @throws IOException
     */
    private void getLogCutadapt(Process proc) throws IOException {

        BufferedReader stdInput = new BufferedReader(new
                InputStreamReader(proc.getInputStream()));
        BufferedReader stdError = new BufferedReader(new
                InputStreamReader(proc.getErrorStream()));

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
     * Method of the class TrimFastq to make some stats of the trimming
     * @param infoTrimPath, a path to store stats in a file
     * @throws IOException
     */
    public void statsLogCutadapt(String infoTrimPath) throws IOException {

        LocalReporter localReporterNumberTimesAdaptor = new LocalReporter();
        BufferedReader infoTrimLeftFile = new BufferedReader(new FileReader(infoTrimPath));
        String line ;
        String oldID="";
        String stackConstructAdaptor="";

        while ((line = infoTrimLeftFile.readLine()) != null) {

            String[] part=line.split("\t");
            String ID = part[0];
            String error = part[1];

            if(error.equals("-1")){

                if(localReporterNumberTimesAdaptor.getCounterNames("Construction").contains("no_adaptor_found")){

                    localReporterNumberTimesAdaptor.incrCounter("Construction", "no_adaptor_found", 1);
                }else{
                    localReporterNumberTimesAdaptor.setCounter("Construction", "no_adaptor_found", 1);
                }
                continue;
            }

            String nameAdaptor=part[7];

            if(!ID.equals(oldID)){

                if(stackConstructAdaptor.equals("")){
                    if(localReporterNumberTimesAdaptor.getCounterNames("Construction").contains(nameAdaptor)){

                        localReporterNumberTimesAdaptor.incrCounter("Construction", nameAdaptor, 1);

                    }else{
                        localReporterNumberTimesAdaptor.setCounter("Construction",nameAdaptor,1);
                    }
                }else {

                    if (stackConstructAdaptor.contains(" + ")) {
                        localReporterNumberTimesAdaptor.incrCounter("Construction", stackConstructAdaptor, 1);
                    }else{
                        localReporterNumberTimesAdaptor.incrCounter("Construction", stackConstructAdaptor, 1);
                    }
                    stackConstructAdaptor = "";
                }

            }else{
                if(stackConstructAdaptor.equals("")){
                    stackConstructAdaptor+=nameAdaptor;
                }else{
                    stackConstructAdaptor+=" + "+nameAdaptor;
                }
            }
            oldID=ID;
        }

        // Analyze counter group Construction

        for(String constructionAdaptor : localReporterNumberTimesAdaptor.getCounterNames("Construction")){
            System.out.println(constructionAdaptor+" :  "+localReporterNumberTimesAdaptor.getCounterValue("Construction", constructionAdaptor));
        }
    }



}
