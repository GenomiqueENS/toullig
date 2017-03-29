package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.toullig.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by birer on 29/03/17.
 */
public class PerfectAlgorithm {

    private HashMap<String, String[]> fastqHash;
    private File nameOutputFastq;
    private int addIndexOutlier;
    private File adaptorFile;
    private int seedMismatchesTrimmomatic;
    private int palindromeClipThresholdTrimmomatic;
    private int simpleClipThreshold;
    private boolean processCutadapt;
    private boolean processTrimmomatic;
    private BufferedWriter fastaFileLeftOutlier;
    private BufferedWriter fastaFileRightOutlier;

    public PerfectAlgorithm(HashMap fastqHash, File nameOutputFastq, int addIndexOutlier,File adaptorFile, int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic, int simpleClipThreshold, boolean processCutadapt, boolean processTrimmomatic, String pathFastaFileLeftOutlier, String pathFastaFileRightOutlier) throws IOException {

        this.fastqHash=fastqHash;
        this.nameOutputFastq=nameOutputFastq;
        this.addIndexOutlier=addIndexOutlier;
        this.adaptorFile=adaptorFile;
        this.seedMismatchesTrimmomatic=seedMismatchesTrimmomatic;
        this.palindromeClipThresholdTrimmomatic=palindromeClipThresholdTrimmomatic;
        this.simpleClipThreshold=simpleClipThreshold;
        this.processCutadapt=processCutadapt;
        this.processTrimmomatic=processTrimmomatic;
        this.fastaFileLeftOutlier= new BufferedWriter(new FileWriter(pathFastaFileLeftOutlier));
        this.fastaFileRightOutlier= new BufferedWriter(new FileWriter(pathFastaFileRightOutlier));


    }

    /**
     * Method of the class TrimModes to trim sequence to create sequences files for cutadapt.
     * @return a HasMap with values update for trimming
     * @throws IOException
     * @throws InterruptedException
     */
    public HashMap trimOutlierWithPMethod() throws IOException, InterruptedException {

        UtilsTrimming utilsTrimming = new UtilsTrimming();

        int count1=0;
        int count2=0;
        int count3=0;
        int count4=0;
        int count5=0;
        int count6=0;
        int count7=0;

        BufferedWriter BufferedWriterFastqOutput = new BufferedWriter(new FileWriter(this.nameOutputFastq));

        Utils utils = new Utils();

        for (String ID : this.fastqHash.keySet()) {

            int lengthOutlierEnd=0;
            int lengthOutlierBegin=0;
            count1++;
            String[] tabValue = this.fastqHash.get(ID);
            String sequence = tabValue[0];
            String score = tabValue[1];
            String cigar = tabValue[2];
            //trim by CIGAR
            if(!cigar.equals("*")){

                //Repere Soft and Hard clipping
                int beginIndexExtremite=cigar.indexOf("S");
                if(beginIndexExtremite>cigar.indexOf("H") && tabValue[2].indexOf("H")>=0){
                    beginIndexExtremite=cigar.indexOf("H");
                }

                if(beginIndexExtremite!=-1 && beginIndexExtremite+1!=cigar.length() && beginIndexExtremite<=6){

                    lengthOutlierBegin = Integer.parseInt(cigar.substring(0,beginIndexExtremite))+this.addIndexOutlier;
                    tabValue[3]=""+lengthOutlierBegin;

                    count2++;
                }else{
                    //System.out.println("CIGAR problem deb: "+tabValue[2]);
                    lengthOutlierBegin=0;
                }

                //Repere Soft and Hard clipping
                int endIndexExtremite=cigar.lastIndexOf("S");
                if(endIndexExtremite<cigar.lastIndexOf("H")){
                    endIndexExtremite=cigar.lastIndexOf("H");
                }
                if(endIndexExtremite==-1){
                    endIndexExtremite=0;
                }

                if(endIndexExtremite+1==cigar.length() && cigar.length()-endIndexExtremite<5 ){

                    lengthOutlierEnd = Integer.parseInt(cigar.substring(cigar.lastIndexOf("M")+1,endIndexExtremite))+this.addIndexOutlier;
                    tabValue[4]=""+lengthOutlierEnd;
                    count3++;
                }else{
                    //System.out.println("CIGAR problem end: "+tabValue[2]);
                    lengthOutlierEnd=0;
                }

                if(this.processCutadapt){
                    utilsTrimming.writeOutliers(lengthOutlierBegin,lengthOutlierEnd,sequence,ID, utils, this.fastaFileLeftOutlier, this.fastaFileRightOutlier);
                }

                if(this.processTrimmomatic){

                    String leftOutlierSequence= utilsTrimming.getOutlierLeftSequence(lengthOutlierBegin,sequence);
                    String rightOutlierSequence= utilsTrimming.getOutlierRightSequence(lengthOutlierEnd,sequence);
                    String leftOutlierScore= utilsTrimming.getOutlierLeftScore(lengthOutlierBegin,score);
                    String rightOutlierScore= utilsTrimming.getOutlierRightScore(lengthOutlierEnd,score);

                    String mainSequence = sequence.substring(lengthOutlierBegin,sequence.length()-lengthOutlierEnd);
                    String mainScore = score.substring(lengthOutlierBegin,score.length()-lengthOutlierEnd);

                    TrimWithTrimmomatic trimmingTrimmomatic= new TrimWithTrimmomatic(this.adaptorFile, this.seedMismatchesTrimmomatic, this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold);

                    String rigthTrimSequence=trimmingTrimmomatic.TrimmomaticTrim(rightOutlierSequence,rightOutlierScore);
                    String leftTrimSequence=trimmingTrimmomatic.TrimmomaticTrim(utils.reverse(leftOutlierSequence),utils.reverse(leftOutlierScore));

                    String rigthTrimScore=score.substring(mainScore.length(),mainScore.length()+rightOutlierSequence.length());
                    String leftTrimScore=score.substring(leftOutlierScore.length()-leftTrimSequence.length(),leftOutlierScore.length());

                    String sequenceTranscript=leftTrimSequence+mainSequence+rigthTrimSequence;
                    String scoreTranscript = leftTrimScore+mainScore+rigthTrimScore;

                    utils.writeFastq( ID, sequenceTranscript, scoreTranscript, BufferedWriterFastqOutput);

                }
            }else{
                count4++;
            }
            if(tabValue[5].equals("16")){
                count5++;
            }
            if(tabValue[5].equals("0")){
                count6++;
            }
            if(tabValue[5].equals("4")){
                count7++;
            }
        }


        this.fastaFileLeftOutlier.close();
        this.fastaFileRightOutlier.close();


        System.out.println("Nombre total de sequence dans le SAM : "+count1);
        System.out.println("Nombre de sequence CIGAR '*': "+count4);
        System.out.println("Nombre de QFlag '16': "+count5);
        System.out.println("Nombre de QFlag '0': "+count6);
        System.out.println("Nombre de QFlag '4': "+count7);
        System.out.println("Nombre d'Outlier Debut trouvé: "+count2);
        System.out.println("Nombre d'Outlier Fin trouvé: "+count3);
        if(this.processTrimmomatic) {
            BufferedWriterFastqOutput.close();
        }
        return this.fastqHash;
    }
}
