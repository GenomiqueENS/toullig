package fr.ens.biologie.genomique.toullig;

import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;
import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.Sequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;


import java.io.*;
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
    private String PCRPrimer;
    private String strandSwitching;
    private BufferedWriter fastaFileLeftOutlier;
    private BufferedWriter fastaFileRightOutlier;
    private Alphabet alphabet = AMBIGUOUS_DNA_ALPHABET;

    public TrimFastq(File samFile, File fastqFile, File nameOutputFastq) throws IOException, InterruptedException {

        String pathOutputTrimLeftFasta= "/home/birer/Bureau/nanoporetools/output/outputFastaFileLeftOutlier.fasta";
        String pathOutputTrimRightFasta= "/home/birer/Bureau/nanoporetools/output/outputFastaFileRightOutlier.fasta";
        String pathFastaFileLeftOutlier = "/home/birer/Bureau/nanoporetools/output/fastaFileLeftOutlier.fasta";
        String pathFastaFileRightOutlier = "/home/birer/Bureau/nanoporetools/output/fastaFileRightOutlier.fasta";


        //Thread.sleep(20000);

        //readAdaptorRTFile(new BufferedReader(new FileReader("/home/birer/Bureau/nanoporetools/config_files/adaptor_RT_sequence.txt")));
        InputStream samInputStream =new FileInputStream(samFile);
        readSamFile(samInputStream);

        readFastqFile(fastqFile);

        this.fastaFileLeftOutlier= new BufferedWriter(new FileWriter(pathFastaFileLeftOutlier));
        this.fastaFileRightOutlier= new BufferedWriter(new FileWriter(pathFastaFileRightOutlier));

        System.out.println("Trim begin !");

        trimOutlier1();

        int lenWindows=15;
        double threshold=0.8;
        //trimOutlier2(lenWindows,threshold);

        this.fastaFileRightOutlier.close();
        this.fastaFileLeftOutlier.close();

        System.out.println("Begin use cutadapt !");

        String adaptorRT="CTTGCCTGTCGCTCTATCTTCTTTTTVN";
        String adaptorSwithStrand="TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG";
        String strandLeft="-g";
        String strandRight="-a";
        double errorRate=0.5;

        cutadaptTrim(pathFastaFileLeftOutlier,pathOutputTrimLeftFasta,adaptorRT,adaptorSwithStrand,strandLeft,errorRate);
        cutadaptTrim(pathFastaFileRightOutlier,pathOutputTrimRightFasta,adaptorRT,adaptorSwithStrand,strandRight,errorRate);

        mergeTrimOutlier(pathOutputTrimLeftFasta,pathOutputTrimRightFasta,nameOutputFastq);
    }


    /**
     * Method of the class TrimFastq to trim sequence to create sequences files for cutadapt.
     * @throws IOException
     * @throws InterruptedException
     */
    public void trimOutlier1() throws IOException, InterruptedException {

        int count1=0;
        int count2=0;
        int count3=0;
        int count4=0;
        int count5=0;
        int count6=0;
        int count7=0;

        for (String key : this.fastqHash.keySet()) {

            int lengthOutlierEnd=0;
            int lengthOutlierBegin=0;
            count1++;
            String[] tabValue = this.fastqHash.get(key);
            String sequence = tabValue[0];
            //trim by CIGAR
            if(!tabValue[2].equals("*")){

                //Repere Soft and Hard clipping
                int beginIndexExtremite=tabValue[2].indexOf("S");
                if(beginIndexExtremite>tabValue[2].indexOf("H") && tabValue[2].indexOf("H")>=0){
                    beginIndexExtremite=tabValue[2].indexOf("H");
                }

                if(beginIndexExtremite!=-1 && beginIndexExtremite+1!=tabValue[2].length() && beginIndexExtremite<=6){
                    tabValue[3]=tabValue[2].substring(0,beginIndexExtremite);
                    lengthOutlierBegin = Integer.parseInt(tabValue[2].substring(0,beginIndexExtremite));

                    count2++;
                }else{
                    //System.out.println("CIGAR problem deb: "+tabValue[2]);
                    lengthOutlierBegin=0;
                }

                //Repere Soft and Hard clipping
                int endIndexExtremite=tabValue[2].lastIndexOf("S");
                if(endIndexExtremite<tabValue[2].lastIndexOf("H")){
                    endIndexExtremite=tabValue[2].lastIndexOf("H");
                }
                if(endIndexExtremite==-1){
                    endIndexExtremite=0;
                }

                if(endIndexExtremite+1==tabValue[2].length() && tabValue[2].length()-endIndexExtremite<5 ){
                    tabValue[4]=tabValue[2].substring(tabValue[2].lastIndexOf("M")+1,endIndexExtremite);
                    lengthOutlierEnd = Integer.parseInt(tabValue[2].substring(tabValue[2].lastIndexOf("M")+1,endIndexExtremite));
                    count3++;
                }else{
                    //System.out.println("CIGAR problem end: "+tabValue[2]);
                    lengthOutlierEnd=0;
                }

                writeOutlier(lengthOutlierBegin,lengthOutlierEnd,sequence,key,this.fastaFileLeftOutlier,this.fastaFileRightOutlier);



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

        System.out.println("Nombre total de sequence mapped: "+count1);
        System.out.println("Nombre de sequence CIGAR '*': "+count4);
        System.out.println("Nombre de QFlag '16': "+count5);
        System.out.println("Nombre de QFlag '0': "+count6);
        System.out.println("Nombre de QFlag '4': "+count7);
        System.out.println("Nombre d'Outlier Debut bien trouvé: "+count2);
        System.out.println("Nombre d'Outlier Fin bien trouvé: "+count3);
    }

    /**
     * Method of the class TrimFastq to trim sequence to create sequences files for cutadapt with the side-windows method.
     * @param lenWindows, the length of the window
     * @param threshold, the threshold
     * @throws IOException
     * @throws InterruptedException
     */
    public void trimOutlier2(int lenWindows,double threshold) throws IOException, InterruptedException {

        Pattern pattern1 = Pattern.compile("(([0-9]*[A-Z]).*)");
        Pattern pattern2 = Pattern.compile("([0-9]*)(.)");
        int longueurSequenceCigar=0;

        for (String key : this.fastqHash.keySet()) {

            String sequenceCigarBinaire="";
            String[] tabValue = this.fastqHash.get(key);
            String sequence = tabValue[0];
            String CIGAR=tabValue[2];

            if(!CIGAR.equals("*")){

                while(CIGAR.length()!=0){

                    Matcher matcher1 = pattern1.matcher(CIGAR);
                    boolean b = matcher1.matches();
                    // si recherche fructueuse
                    if(b) {

                        // pour chaque groupe
                        String oneCigare=matcher1.group(2);
                        CIGAR=CIGAR.substring(matcher1.group(2).length(),CIGAR.length());
                        Matcher matcher2 = pattern2.matcher(oneCigare);
                        boolean b2 = matcher2.matches();
                        // si recherche fructueuse
                        if(b2){

                            longueurSequenceCigar=Integer.parseInt(matcher2.group(1))-1;
                            for(int i=0;i<=longueurSequenceCigar;i++){


                                if(matcher2.group(2).equals("M")){
                                    sequenceCigarBinaire=sequenceCigarBinaire+1;
                                    continue;
                                }
                                if(matcher2.group(2).equals("N") || matcher2.group(2).equals("D")){
                                    continue;
                                }
                                if(!matcher2.group(2).equals("N") || !matcher2.group(2).equals("M")){
                                    sequenceCigarBinaire=sequenceCigarBinaire+0;
                                }

                            }
                        }
                    }
                }

                int lengthOutlierBegin=sideWindowsLeft(sequenceCigarBinaire,lenWindows,threshold);
                int lengthOutlierEnd=sideWindowsRight(sequenceCigarBinaire,lenWindows,threshold);
                tabValue[3]=""+lengthOutlierBegin;
                tabValue[4]=""+lengthOutlierEnd;

                writeOutlier(lengthOutlierBegin,lengthOutlierEnd,sequence,key,this.fastaFileLeftOutlier,this.fastaFileRightOutlier);
            }
        }
    }

    //
    //Merge outlier
    //

    /**
     * Method of the class TrimFastq to merge the results of cutadapt with the trim sequence.
     * @param pathOutputTrimLeftFasta, the path of the left output
     * @param pathOutputTrimRightFasta, the path of the right output
     * @throws IOException
     */
    public void mergeTrimOutlier(String pathOutputTrimLeftFasta,String pathOutputTrimRightFasta, File nameOutputFastq) throws IOException {

        BufferedWriter fastqTrimBufferedWritter= new BufferedWriter(new FileWriter(nameOutputFastq));
        HashMap<String, String> fastaLeftHash = new HashMap<String, String>();
        HashMap<String, String> fastaRightHash = new HashMap<String, String>();

        //Read of the left outlier fasta output by cutadapt
        try{
            BufferedReader outputTrimLeftFastaReader = new BufferedReader(new FileReader(pathOutputTrimLeftFasta));
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
            BufferedReader outputTrimRightFastaFile = new BufferedReader(new FileReader(pathOutputTrimRightFasta));
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
            String leftSequence="";
            String rightSequence="";

            int lengthBeginOutlier=0;
            int lengthEndOutlier=0;

            String[] tabValue = this.fastqHash.get(ID);


            if(tabValue[3].equals("")){
                lengthBeginOutlier=0;
            }else{
                lengthBeginOutlier=Integer.parseInt(tabValue[3]);
            }
            if(tabValue[4].equals("")){
                lengthEndOutlier=0;
            }else{
                lengthEndOutlier=Integer.parseInt(tabValue[4]);
            }

            String mainSequenceWithoutOutlier=tabValue[0].substring(lengthBeginOutlier,tabValue[0].length()-lengthEndOutlier);

            if(fastaLeftHash.get(ID)==null){
                leftSequence="";
            }else{
                leftSequence=fastaLeftHash.get(ID);
            }
            if(fastaLeftHash.get(ID)==null){
                rightSequence="";
            }else{
                rightSequence=fastaRightHash.get(ID);
            }

            String scoreTrim=tabValue[1].substring(lengthBeginOutlier-leftSequence.length(),(tabValue[1].length()-lengthEndOutlier)+rightSequence.length());
            String sequenceTrim =leftSequence+mainSequenceWithoutOutlier+rightSequence;

            writeFastq(ID,sequenceTrim,scoreTrim,fastqTrimBufferedWritter);

        }
        fastqTrimBufferedWritter.close();
    }

    //
    //Cutadapt
    //

    /**
     * Method of the class TrimFastq to execute cutadapt.
     * @param pathFastaFileOutlier, the path of the fasta file with outlier to trim
     * @param pathOutputTrimFasta, the path of the trimmed reads in fasta file
     * @param adaptorRT, the sequence of the adaptor RT
     * @param adaptorSwithStrand, the sequence of the adaptor SwitchStrand
     * @param strand, the cutadapt strand of the adaptor
     * @throws IOException
     * @throws InterruptedException
     */
    public void cutadaptTrim(String pathFastaFileOutlier,String pathOutputTrimFasta,String adaptorRT,String adaptorSwithStrand,String strand, double errorRate) throws IOException, InterruptedException {

        try {
            String reverseComplementAdaptorRT=strand+" reverse_complement_RT_adaptor="+reverseComplement(adaptorRT,this.alphabet);
            String reverseComplementAdaptorSwithStrand=strand+" reverse_complement_Switch_Strand_RT_adaptor="+reverseComplement(adaptorSwithStrand,this.alphabet);

            String complementAdaptorRT=strand+" complement_RT_adaptor="+complement(adaptorRT,this.alphabet);
            String complementAdaptorSwithStrand=strand+" complement_Switch_Strand_RT_adaptor="+complement(adaptorSwithStrand,this.alphabet);

            StringBuffer reverseAdaptorRTStringBuffer = new StringBuffer(adaptorRT);
            StringBuffer reverseAdaptorSwithStrandStringBuffer = new StringBuffer(adaptorSwithStrand);
            String reverseAdaptorRT= strand+" reverse_RT_adaptor="+reverseAdaptorRTStringBuffer.reverse().toString();
            String reverseAdaptorSwithStrand= strand+" reverse_Switch_Strand_RT_adaptor="+reverseAdaptorSwithStrandStringBuffer.reverse().toString();

            adaptorRT=strand+" RT_adaptor="+adaptorRT;
            adaptorSwithStrand=strand+" Switch_Strand_RT_adaptor="+adaptorSwithStrand;

            //--quiet
            ProcessBuilder pb = new ProcessBuilder("/bin/bash", "-c", "cutadapt "+adaptorRT+" "+adaptorSwithStrand+" "+reverseAdaptorRT+" "+reverseAdaptorSwithStrand+" "+complementAdaptorRT+" "+complementAdaptorSwithStrand+" "+reverseComplementAdaptorRT+" "+reverseComplementAdaptorSwithStrand+""+" --error-rate="+errorRate+" --overlap=6 --times=8 --match-read-wildcards --format=fasta "+pathFastaFileOutlier+" > "+pathOutputTrimFasta);

            //System.out.println(pb.command());
            pb.redirectErrorStream(true);
            Process proc = pb.start(); // Start the process.

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

            //int exitproc=proc.waitFor();
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

    //
    // get Outlier
    //

    /**
     * Method of the class TrimFastq to write the outliers (3' and 5') in fasta files.
     * @param lengthOutlierBegin, the length of the outlier 3'
     * @param lengthOutlierEnd, the length of the outlier 5'
     * @param sequence, the sequence of the read
     * @param ID, the ID of the read
     * @param fastaFileLeftOutlier, the file (.fasta) of the outlier 3'
     * @param fastaFileRightOutlier, the file (.fasta) of the outlier 5'
     * @throws IOException
     */
    public void writeOutlier(int lengthOutlierBegin, int lengthOutlierEnd, String sequence, String ID, BufferedWriter fastaFileLeftOutlier, BufferedWriter fastaFileRightOutlier) throws IOException {

        String leftOutlierSequence=getOutlierLeft(lengthOutlierBegin,sequence);
        String rightOutlierSequence= getOutlierRight(lengthOutlierEnd,sequence);
        writeFasta(leftOutlierSequence,ID,fastaFileLeftOutlier);
        writeFasta(rightOutlierSequence,ID,fastaFileRightOutlier);
    }

    /**
     * Method of the class TrimFastq to obtain the sequence of the 3' (Left) outlier.
     * @param lengthOutlierBegin, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the left outlier
     */
    public String getOutlierLeft(int lengthOutlierBegin, String sequence){

        return sequence.substring(0,lengthOutlierBegin);
    }

    /**
     * Method of the class TrimFastq to obtain the sequence of the 5' (Right) outlier.
     * @param lengthOutlierEnd, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the right outlier
     */
    public String getOutlierRight(int lengthOutlierEnd, String sequence){

        return sequence.substring(sequence.length()-lengthOutlierEnd,sequence.length());
    }

    //
    // get index by side window
    //

    /**
     * Method of the class TrimFastq to the length of the left outlier with a binair CIGAR sequence.
     * @param sequenceCigarBinaire, a binaire CIGAR sequence
     * @param lenWindows, the length of the window
     * @param threshold, the threshold to separate the outlier to the main sequence
     * @return int, the length of the outlier
     */
    public int sideWindowsLeft(String sequenceCigarBinaire,int lenWindows,double threshold){

        if(sequenceCigarBinaire==null){
            return 0;
        }
        String windows="";
        int length=0;
        for(int i = 0;i<=sequenceCigarBinaire.length();i++){
            if(i==(sequenceCigarBinaire.length()-lenWindows-2)){
                break;
            }
            windows=sequenceCigarBinaire.substring(i,i+lenWindows);
            if(sumWindowCIGAR(windows)>=threshold){
                length=i+lenWindows;
                break;
            }
        }
        return length;
    }

    /**
     * Method of the class TrimFastq to the length of the right outlier with a binair CIGAR sequence.
     * @param sequenceCigarBinaire, a binaire CIGAR sequence
     * @param lenWindows, the length of the window
     * @param threshold, the threshold to separate the outlier to the main sequence
     * @return int, the length of the outlier
     */
    public int sideWindowsRight(String sequenceCigarBinaire,int lenWindows,double threshold){

        if(sequenceCigarBinaire==null){
            return 0;
        }
        String windows="";
        int length=0;
        for(int i = sequenceCigarBinaire.length();i>=0;i--){
            if(i==lenWindows){
                break;
            }

            windows=sequenceCigarBinaire.substring(i-lenWindows,i);
            if(sumWindowCIGAR(windows)>=threshold){
                if(i>=sequenceCigarBinaire.length()-lenWindows){
                    length=sequenceCigarBinaire.length()-i;
                }else{
                    length=sequenceCigarBinaire.length()-i-lenWindows;
                }
                break;
            }
        }
        return length;
    }


    /**
     * Method of the class TrimFastq to sumWindowCIGAR a sequence CIGAR encode in 1/0.
     * @param windows, a String of CIGAR sequence encode in 1/0
     * @return double, the sumWindowCIGAR of the CIGAR window
     */
    public double sumWindowCIGAR(String windows){
        String[] arrayWindows = windows.split("");
        int sum=0;
        for(int i = 0;i<=arrayWindows.length-1;i++){
            if(arrayWindows[i].equals("1")){
                sum+=1;
            }
        }
        return (double)sum/windows.length();
    }

    //
    // Write
    //

    /**
     * Method of the class TrimFastq to write a sequence to the fasta format.
     * @param sequence, the sequence of the read
     * @param ID, the ID of the read
     * @param fastaFile, the file to write fasta sequence
     * @throws IOException
     */
    public void writeFasta(String sequence, String ID, BufferedWriter fastaFile) throws IOException {

        try {
            fastaFile.write(">"+ID+"\n");
            for (int i = 0; i <= sequence.length(); i = i + 60) {
                if (i + 60 >= sequence.length()) {
                    fastaFile.write(sequence.substring(i, sequence.length()) + "\n");
                    break;
                }
                fastaFile.write(sequence.substring(i, i + 60) + "\n");
            }
        }catch (IOException e){
            e.printStackTrace();
        }
    }


    /**
     * Method of the class TrimFastq to write a sequence to the fastq format.
     * @param ID, the ID of the read
     * @param sequenceTrim, the sequence of the read
     * @param scoreTrim, the score of the read
     * @param fastqTrimBufferedWritter, the file to write fastq sequence
     * @throws IOException
     */
    public void writeFastq(String ID,String sequenceTrim,String scoreTrim, BufferedWriter fastqTrimBufferedWritter)throws IOException{

        fastqTrimBufferedWritter.write("@"+ID+"\n");
        fastqTrimBufferedWritter.write(sequenceTrim+"\n");
        fastqTrimBufferedWritter.write("+\n");
        fastqTrimBufferedWritter.write(scoreTrim+"\n");
    }

    //
    // Read Files
    //

    /**
     * Method of the class TrimFastq to read and get the adaptor in the adaptor file.
     * @param adaptorFile, the file who contains adaptor
     * @throws IOException
     */
    public void readAdaptorRTFile(BufferedReader adaptorFile)throws IOException{

        int i =0;
        String line;
        while ((line = adaptorFile.readLine()) != null) {
            if(i==1){
                this.PCRPrimer=line;
            }
            if(i==3){
                this.strandSwitching=line;
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
    public void readSamFile(InputStream samInputStream) throws IOException {

        try(final SamReader inputSam =
                    SamReaderFactory.makeDefault().open(SamInputResource.of(samInputStream))){

            String lengthBeginOutlier="";
            String lengthEndOutlier="";
            String sequence="";
            String score="";
            for (SAMRecord samRecord : inputSam) {

                String CIGAR=samRecord.getCigarString();
                String QFlag = ""+samRecord.getFlags();
                String ID =samRecord.getReadName();

                if(this.fastqHash.containsKey(ID)) {
                    this.fastqHash.put(ID, new String[]{sequence, score, CIGAR, lengthBeginOutlier, lengthEndOutlier, QFlag});
                }else {
                    this.fastqHash.put(ID, new String[]{sequence, score, CIGAR, lengthBeginOutlier, lengthEndOutlier, QFlag});
                }
            }
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
    public void readFastqFile(File fastqBufferedReader)throws IOException{

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
    //  utils
    //


    /**
     * Get the sequence as the complement. This method work only with
     * A,T,G and C bases.
     * @param sequence sequence to reverse complement
     * @param alphabet alphabet of the sequence to reverse complement
     * @return the reverse complement sequence
     */
    public static final String complement(final String sequence,
                                          final Alphabet alphabet) {

        if (sequence == null || alphabet == null) {
            return null;
        }

        String s =
                reverseComplement(sequence,alphabet);
        final char[] array = s.toCharArray();
        final int len = array.length;
        final StringBuilder sb = new StringBuilder(len);

        for (int i = len-1; i >=  0; i--) {
            sb.append(array[i]);
        }
        return sb.toString();
    }
}

