/**
 * 
 */
package com.sunset;

import java.net.URL;
import java.util.Vector;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.apache.axis.client.Call;
import org.apache.axis.client.Service;
import org.apache.axis.message.SOAPBodyElement;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * VCCLAB web service provides calculation of logP and logS this class is a web service client
 * @author Oleg Ursu
 *
 */
public class VCCLab {

	private String smiles;
	private double logP;
	private double logS;
	
	public VCCLab() {
		
	};
	
	public void setSMILES(String smiles) {
		this.smiles = smiles;
	}
	
	@SuppressWarnings("rawtypes")
	public void run() throws Exception {
		Service service = new Service();
		Call call = (Call)service.createCall();
		String endpoint = "http://www.vcclab.org/web/services/ALOGPS";
		call.setTargetEndpointAddress(new URL(endpoint));
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document doc = builder.newDocument();
		Element rootElement = doc.createElement("MOLECULES");
		rootElement.setAttribute("FORMAT", "SMILES");
		Element molecule = doc.createElement("MOLECULE");
		molecule.appendChild(doc.createTextNode(smiles));
		rootElement.appendChild(molecule);
		doc.appendChild(rootElement);
		SOAPBodyElement[] input = new SOAPBodyElement[1];
		input[0] = new SOAPBodyElement(doc.getDocumentElement());
		Vector elems = (Vector) call.invoke(input);
	    SOAPBodyElement element = (SOAPBodyElement) elems.get(0);
	    Element e = element.getAsDOM();
	    Element molElement = (Element) e.getElementsByTagName("MOLECULE").item(0);
	    logP = Double.parseDouble(((Element)molElement.getElementsByTagName("LOGP").item(0)).getFirstChild().getNodeValue());
	    logS = Double.parseDouble(((Element)molElement.getElementsByTagName("LOGS").item(0)).getFirstChild().getNodeValue());
	}
	
	public double getLogP() {
		return logP;
	}
	
	public double getLogS() {
		return logS;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		VCCLab vcc = new VCCLab();
		vcc.setSMILES("CCO");
		vcc.run();
		System.out.println("logP = " + vcc.getLogP());
		System.out.println("logS = " + vcc.getLogS());
	}

}
