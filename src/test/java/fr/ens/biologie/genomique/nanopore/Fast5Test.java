package fr.ens.biologie.genomique.nanopore;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Paths;
import java.util.Date;

import org.junit.Test;

import fr.ens.biologie.genomique.nanopore.Fast5.ChemistryVersion;
import fr.ens.biologie.genomique.nanopore.Fast5.Status;
import fr.ens.biologie.genomique.nanopore.Fast5.Type;
import fr.ens.biologie.genomique.nanopore.Fast5.Version;



public class Fast5Test {
	String file1 = "/alexander_PC_20161027_R9-4_1D.fast5";
	String file2 = "/dnacpc14_20160617_R7_2D_prebasecalling.fast5";
	String file3 = "/dnacpc14_20160617_R7_2D.fast5";
	String file4 = "/dnacpc14_20161011_R9_2D_prebasecalling.fast5";
	String file5 = "/dnacpc14_20161011_R9_2D.fast5";

	/**
	 * Read an input stream.
	 * 
	 * @throws MalformedURLException
	 * @throws IOException
	 */
	private static String readInputStream(String path) throws IOException {

		InputStream is = Fast5Test.class.getResourceAsStream(path);
		StringBuilder sb = new StringBuilder();
		try (BufferedReader br = new BufferedReader(new InputStreamReader(is))) {
			String s = null;
			while ((s = br.readLine()) != null) {
				sb.append(s);
				sb.append('\n');
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sb.toString();
	}
	
	private File getResourceAsFile(String path){
		URL resource = Fast5Test.class.getResource(path);
		try {
			return Paths.get(resource.toURI()).toFile();
		} catch (URISyntaxException e) {
			return null;
		}
	}
	
	@Test
	public void testReadVersion() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(Version.V1_1, testf1.getVersion());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(null, testf2.getVersion());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(Version.V1_1, testf3.getVersion());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(null, testf4.getVersion());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(Version.V1_1, testf5.getVersion());
		testf5.close();

	}

	@Test
	public void testReadStatus() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(Status.AFTER_BASECALLING, testf1.getStatus());
		assertNotEquals(Status.PRE_BASECALLING, testf1.getStatus());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(Status.PRE_BASECALLING, testf2.getStatus());
		assertNotEquals(Status.AFTER_BASECALLING, testf2.getStatus());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(Status.AFTER_BASECALLING, testf3.getStatus());
		assertNotEquals(Status.PRE_BASECALLING, testf3.getStatus());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(Status.PRE_BASECALLING, testf4.getStatus());
		assertNotEquals(Status.AFTER_BASECALLING, testf4.getStatus());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(Status.AFTER_BASECALLING, testf5.getStatus());
		assertNotEquals(Status.PRE_BASECALLING, testf5.getStatus());
		testf5.close();

	}

	@Test
	public void testReadType() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(Type.TYPE_1D, testf1.getType());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(null, testf2.getType());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(Type.TYPE_2D, testf3.getType());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(null, testf4.getType());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(Type.TYPE_2D, testf5.getType());
		testf5.close();

	}

	@Test
	public void testReadRVersion() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(ChemistryVersion.R9_4, testf1.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R9, testf1.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R7, testf1.getChemistryVersion());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(ChemistryVersion.R7, testf2.getChemistryVersion());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(ChemistryVersion.R7, testf3.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R9, testf3.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R9_4, testf3.getChemistryVersion());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(ChemistryVersion.R9, testf4.getChemistryVersion());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(ChemistryVersion.R9, testf5.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R9_4, testf5.getChemistryVersion());
		assertNotEquals(ChemistryVersion.R7, testf5.getChemistryVersion());
		testf5.close();

	}

	//
	//
	// macro
	//
	//

	@Test
	public void testIsBarcoded() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertFalse(testf1.isBarcoded());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertFalse(testf2.isBarcoded());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertFalse(testf3.isBarcoded());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertFalse(testf4.isBarcoded());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertFalse(testf5.isBarcoded());
		testf5.close();
	}

	@Test
	public void testIsBasecalled() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(true, testf1.isBasecalled());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(false, testf2.isBasecalled());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(true, testf3.isBasecalled());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(false, testf4.isBasecalled());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(true, testf5.isBasecalled());
		testf5.close();
	}

	@Test
	public void testIs2D() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(false, testf1.is2D());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(false, testf2.is2D());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(true, testf3.is2D());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(false, testf4.is2D());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(true, testf5.is2D());
		testf5.close();
	}

	//
	//
	// primary information getters
	//
	//

	//
	// Raw Group Information
	//

	// @Test
	// public void testGetElectricalSignal() {
	// Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	// assertEquals(998, testf1.getElectricalSignal()[0]);
	// testf1.close();
	//
	// Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	// assertEquals(51.24900931222098, testf2.getElectricalSignal()[0]);
	// testf2.close();
	//
	// Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	// assertEquals(14614880, testf3.getElectricalSignal()[0]);
	// testf3.close();
	//
	// Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	// assertEquals(1390, testf4.getElectricalSignal()[0]);
	// testf4.close();
	//
	// Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	// assertEquals(1219, testf5.getElectricalSignal()[0]);
	// testf5.close();
	// }

	//
	// UniqueGlobalkey Group Information
	//
	@Test
	public void testGetNumMinION() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("MN17351", testf1.getNumMinION());
		assertNotEquals(null, testf1.getNumMinION());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("MN16592", testf2.getNumMinION());
		assertNotEquals(null, testf2.getNumMinION());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("MN16592", testf3.getNumMinION());
		assertNotEquals(null, testf3.getNumMinION());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("MN17734", testf4.getNumMinION());
		assertNotEquals(null, testf4.getNumMinION());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("MN17734", testf5.getNumMinION());
		assertNotEquals(null, testf5.getNumMinION());
		testf5.close();
	}

	@Test
	public void testGetFlowcellId() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("FAB44188", testf1.getFlowcellId());
		assertNotEquals(null, testf1.getFlowcellId());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("FAA105486", testf2.getFlowcellId());
		assertNotEquals(null, testf2.getFlowcellId());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("FAA105486", testf3.getFlowcellId());
		assertNotEquals(null, testf3.getFlowcellId());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("FAD22491", testf4.getFlowcellId());
		assertNotEquals(null, testf4.getFlowcellId());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("FAD22491", testf5.getFlowcellId());
		assertNotEquals(null, testf5.getFlowcellId());
		testf5.close();
	}

	@Test
	public void testGetFlowcellVersion() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("1.1.15", testf1.getFlowcellVersion());
		assertNotEquals(null, testf1.getFlowcellVersion());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("0.51.3.40 b201605171140", testf2.getFlowcellVersion());
		assertNotEquals(null, testf2.getFlowcellVersion());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("0.51.3.40 b201605171140", testf3.getFlowcellVersion());
		assertNotEquals(null, testf3.getFlowcellVersion());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("1.0.8", testf4.getFlowcellVersion());
		assertNotEquals(null, testf4.getFlowcellVersion());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("1.0.8", testf5.getFlowcellVersion());
		assertNotEquals(null, testf5.getFlowcellVersion());
		testf5.close();
	}

	@Test
	public void testGetdDateExp() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		long Dateint1 = 1477580791;
		Date datef1 = new Date(Dateint1 * 1000);
		assertEquals(datef1, testf1.getdDateExp());
		assertNotEquals(null, testf1.getdDateExp());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		long Dateint2 = 1466159977;
		Date datef2 = new Date(Dateint2 * 1000);
		assertEquals(datef2, testf2.getdDateExp());
		assertNotEquals(null, testf2.getdDateExp());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		long Dateint3 = 1466159977;
		Date datef3 = new Date(Dateint3 * 1000);
		assertEquals(datef3, testf3.getdDateExp());
		assertNotEquals(null, testf3.getdDateExp());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		long Dateint4 = 1476200754;
		Date datef4 = new Date(Dateint4 * 1000);
		assertEquals(datef4, testf4.getdDateExp());
		assertNotEquals(null, testf4.getdDateExp());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		long Dateint5 = 1476200754;
		Date datef5 = new Date(Dateint5 * 1000);
		assertEquals(datef5, testf5.getdDateExp());
		assertNotEquals(null, testf5.getdDateExp());
		testf5.close();
	}

	@Test
	public void testGetProtocolRunId() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("4bb23786-e23f-4d30-9b17-c7ee511fa306", testf1.getProtocolRunId());
		assertNotEquals(null, testf1.getProtocolRunId());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("0cc6bdd9-6487-48be-a9cc-ccad13f25bcf", testf2.getProtocolRunId());
		assertNotEquals(null, testf2.getProtocolRunId());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("0cc6bdd9-6487-48be-a9cc-ccad13f25bcf", testf3.getProtocolRunId());
		assertNotEquals(null, testf3.getProtocolRunId());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("a5dddcb4-2268-4aa1-b44b-641ac2977f1b", testf4.getProtocolRunId());
		assertNotEquals(null, testf4.getProtocolRunId());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("a5dddcb4-2268-4aa1-b44b-641ac2977f1b", testf5.getProtocolRunId());
		assertNotEquals(null, testf5.getProtocolRunId());
		testf5.close();
	}

	@Test
	public void testGetHostname() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("alexander-PC", testf1.getHostname());
		assertNotEquals(null, testf1.getHostname());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("dnacpc14", testf2.getHostname());
		assertNotEquals(null, testf2.getHostname());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("dnacpc14", testf3.getHostname());
		assertNotEquals(null, testf3.getHostname());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("dnacpc14", testf4.getHostname());
		assertNotEquals(null, testf4.getHostname());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("dnacpc14", testf5.getHostname());
		assertNotEquals(null, testf5.getHostname());
		testf5.close();
	}

	@Test
	public void testGetOS() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("Windows 6.1", testf1.getOS());
		assertNotEquals(null, testf1.getOS());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("Windows 6.1", testf2.getOS());
		assertNotEquals(null, testf2.getOS());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("Windows 6.1", testf3.getOS());
		assertNotEquals(null, testf3.getOS());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("Windows 6.1", testf4.getOS());
		assertNotEquals(null, testf4.getOS());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("Windows 6.1", testf5.getOS());
		assertNotEquals(null, testf5.getOS());
		testf5.close();
	}

	@Test
	public void testGetExperimentKit() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("genomic_dna", testf1.getExperimentKit());
		assertNotEquals(null, testf1.getExperimentKit());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals("lambda_burn_in", testf2.getExperimentKit());
		assertNotEquals(null, testf2.getExperimentKit());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("lambda_burn_in", testf3.getExperimentKit());
		assertNotEquals(null, testf3.getExperimentKit());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals("genomic_dna", testf4.getExperimentKit());
		assertNotEquals(null, testf4.getExperimentKit());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("genomic_dna", testf5.getExperimentKit());
		assertNotEquals(null, testf5.getExperimentKit());
		testf5.close();
	}

	@Test
	public void testGetExperimentType() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("customer_qc", testf1.getExperimentType());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(null, testf2.getExperimentType());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(null, testf3.getExperimentType());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(null, testf4.getExperimentType());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("customer_qc", testf5.getExperimentType());
		testf5.close();
	}

	@Test
	public void testGetSampleFrequency() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(4000, testf1.getSampleFrequency());
		assertNotEquals(null, testf1.getSampleFrequency());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(3012, testf2.getSampleFrequency());
		assertNotEquals(null, testf2.getSampleFrequency());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(3012, testf3.getSampleFrequency());
		assertNotEquals(null, testf3.getSampleFrequency());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(4000, testf4.getSampleFrequency());
		assertNotEquals(null, testf4.getSampleFrequency());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(4000, testf5.getSampleFrequency());
		assertNotEquals(null, testf5.getSampleFrequency());
		testf5.close();
	}
	
	@Test
	public void testGetChannelNumber() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		int test = 96;
		assertEquals(test, testf1.getChannelNumber());
		assertNotEquals(null, testf1.getChannelNumber());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(96, testf2.getChannelNumber());
		assertNotEquals(null, testf2.getChannelNumber());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(96, testf3.getChannelNumber());
		assertNotEquals(null, testf3.getChannelNumber());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(96, testf4.getChannelNumber());
		assertNotEquals(null, testf4.getChannelNumber());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(96, testf5.getChannelNumber());
		assertNotEquals(null, testf5.getChannelNumber());
		testf5.close();
	}


	//
	//
	// basecalling information getters
	//
	//
	@Test
	public void testGetNumberRead() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(1000, testf1.getNumberRead());
		assertNotEquals(null, testf1.getNumberRead());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(563, testf2.getNumberRead());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(293, testf3.getNumberRead());
		assertNotEquals(null, testf3.getNumberRead());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(77, testf4.getNumberRead());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(14, testf5.getNumberRead());
		assertNotEquals(null, testf5.getNumberRead());
		testf5.close();
	}

	@Test
	public void testGetBaseCaller() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals("chimaera v1.23.3 | dragonet v1.23.0", testf1.getBaseCaller());
		assertNotEquals(null, testf1.getBaseCaller());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(null, testf2.getBaseCaller());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals("chimaera v1.22.6 | dragonet v1.22.2", testf3.getBaseCaller());
		assertNotEquals(null, testf3.getBaseCaller());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(null, testf4.getBaseCaller());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals("chimaera v1.22.10 | dragonet v1.22.4", testf5.getBaseCaller());
		assertNotEquals(null, testf5.getBaseCaller());
		testf5.close();
	}

	// @Test
	// public void testGetAlignment() {
	// Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	// assertEquals(null, testf1.getAlignment());
	// testf1.close();
	//
	// Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	// assertEquals(null, testf2.getAlignment());
	// testf2.close();
	//
	// Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	// assertEquals(5, testf3.getAlignment());
	// testf3.close();
	//
	// Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	// assertEquals(3, testf4.getAlignment());
	// testf4.close();
	//
	// Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	// assertEquals(null, testf5.getAlignment());
	// assertNotEquals("un alignement ....", testf5.getAlignment());
	// testf5.close();
	// }

	@Test
	public void testGetTemplateLength() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(1193, testf1.getTemplateLength());
		assertNotEquals(null, testf1.getTemplateLength());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(0, testf2.getTemplateLength());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(294, testf3.getTemplateLength());
		assertNotEquals(null, testf3.getTemplateLength());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertNotEquals(null, testf4.getTemplateLength());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(781, testf5.getTemplateLength());
		assertNotEquals(null, testf5.getTemplateLength());
		testf5.close();
	}

	@Test
	public void getComplementeLength() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(0, testf1.getComplementeLength());
		assertNotEquals("1193000000", testf1.getComplementeLength());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(0, testf2.getComplementeLength());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(239, testf3.getComplementeLength());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(0, testf4.getComplementeLength());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(721, testf5.getComplementeLength());
		testf5.close();
	}

	@Test
	public void testGetNumBarcode() {
		Fast5 testf1 = new Fast5(getResourceAsFile(file1));
		assertEquals(null, testf1.getNumBarcode());
		assertNotEquals("BCO5", testf1.getNumBarcode());
		testf1.close();

		Fast5 testf2 = new Fast5(getResourceAsFile(file2));
		assertEquals(null, testf2.getNumBarcode());
		assertNotEquals("BCO5", testf2.getNumBarcode());
		testf2.close();

		Fast5 testf3 = new Fast5(getResourceAsFile(file3));
		assertEquals(null, testf3.getNumBarcode());
		assertNotEquals("BCO5", testf3.getNumBarcode());
		testf3.close();

		Fast5 testf4 = new Fast5(getResourceAsFile(file4));
		assertEquals(null, testf4.getNumBarcode());
		assertNotEquals("BCO5", testf4.getNumBarcode());
		testf4.close();

		Fast5 testf5 = new Fast5(getResourceAsFile(file5));
		assertEquals(null, testf5.getNumBarcode());
		assertNotEquals("BCO5", testf5.getNumBarcode());
		testf5.close();
	}

	//
	//
	// FASTQ getters
	//
	//

	 @Test
	 public void testGetTemplateFastq() throws MalformedURLException, IOException {
	 Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	 String sequence1=readInputStream("/alexander_PC_20161027_R9-4_1D_template.fastq");
	 assertEquals(sequence1, testf1.getTemplateFastq());
	 testf1.close();
	
	 Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	 assertEquals(null, testf2.getTemplateFastq());
	 testf2.close();
	
	 Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	 String sequence3=readInputStream("/dnacpc14_20160617_R7_2D_template.fastq");
	 assertEquals(sequence3, testf3.getTemplateFastq());
	 testf3.close();
	
	 Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	 assertEquals(null, testf4.getTemplateFastq());
	 testf4.close();
	
	 Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	 String sequence5=readInputStream("/dnacpc14_20161011_R9_2D_template.fastq");
	 assertEquals(sequence5, testf5.getTemplateFastq());
	 testf5.close();
	 }
	
	 @Test
	 public void testGetComplementFastq() throws MalformedURLException, IOException {
	 Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	 assertEquals(null, testf1.getComplementFastq());
	 testf1.close();
	
	 Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	 assertEquals(null, testf2.getComplementFastq());
	 testf2.close();
	
	 Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	 String sequence3=readInputStream("/dnacpc14_20160617_R7_2D_complement.fastq");
	 assertEquals(sequence3, testf3.getComplementFastq());
	 testf3.close();
	
	 Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	 assertEquals(null, testf4.getComplementFastq());
	 testf4.close();
	
	 Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	 String sequence5=readInputStream("/dnacpc14_20161011_R9_2D_complement.fastq");
	 assertEquals(sequence5, testf5.getComplementFastq());
	 testf5.close();
	 }
	
	 @Test
	 public void getHairpinFastq() throws MalformedURLException, IOException {
	 Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	 assertEquals(null, testf1.getHairpinFastq());
	 testf1.close();
	
	 Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	 assertEquals(null, testf2.getHairpinFastq());
	 testf2.close();
	
	 Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	 String sequence3=readInputStream("/dnacpc14_20160617_R7_2D_hairpin.fastq");
	 assertEquals(sequence3, testf3.getHairpinFastq());
	 testf3.close();
	
	 Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	 assertEquals(null, testf4.getHairpinFastq());
	 testf4.close();
	
	 Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	 String sequence5=readInputStream("/dnacpc14_20161011_R9_2D_hairpin.fastq");
	 assertEquals(sequence5, testf5.getHairpinFastq());
	 testf5.close();
	 }
	
	 @Test
	 public void testGetBarcodingFastq() {
	 Fast5 testf1 = new Fast5(getResourceAsFile(file1));
	 assertEquals(null, testf1.getBarcodingFastq());
	 testf1.close();
	
	 Fast5 testf2 = new Fast5(getResourceAsFile(file2));
	 assertEquals(null, testf2.getBarcodingFastq());
	 testf2.close();
	
	 Fast5 testf3 = new Fast5(getResourceAsFile(file3));
	 assertEquals(null, testf3.getBarcodingFastq());
	 testf3.close();
	
	 Fast5 testf4 = new Fast5(getResourceAsFile(file4));
	 assertEquals(null, testf4.getBarcodingFastq());
	 testf4.close();
	
	 Fast5 testf5 = new Fast5(getResourceAsFile(file5));
	 assertEquals(null, testf5.getBarcodingFastq());
	 testf5.close();
	 }

}
