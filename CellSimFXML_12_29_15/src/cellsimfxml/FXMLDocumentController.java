/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cellsimfxml;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ResourceBundle;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Slider;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.AnchorPane;
import javafx.scene.paint.Color;
import javafx.scene.text.Text;


/**
 *
 * @author varderes
 */



public class FXMLDocumentController implements Initializable {
   
    @FXML
    private Button quitButton;
    @FXML
    private Button run;
    @FXML
    private Button clear;
    
    @FXML
    private ProgressBar progressBar;
    @FXML
    private Text progressText;
    
    @FXML
    private Text status;
    
    int nx, ny, nz;
    int thresh;
    Task copyWorker;
    
    @FXML
    private Slider sliderX;
    @FXML
    private Slider sliderY;
    @FXML
    private Slider sliderZ;
    @FXML
    private Slider sliderT;
    @FXML
    private Slider sliderSimTime;
    @FXML
    private Slider sliderSR;
    @FXML 
    private Label xLabel;
    @FXML
    private Label yLabel;
    @FXML
    private Label zLabel;
    @FXML
    private Label threshLabel;
    @FXML
    private Label simtimeLabel;
    @FXML
    private Label srloadLabel;

    @FXML 
    private ChoiceBox selection;
    
    ObservableList<String> selections = FXCollections.observableArrayList("Run", "Run & Save Plots", "Run w/ Mathematica");
    
    @FXML 
    private ImageView heartViewer;
    Image heartbeating = new Image(getClass().getResourceAsStream("heartbeating.gif"));
       
    
    @FXML
    private Text warning1;
    @FXML
    private Text optionWarning;
    
    
    AnchorPane AnchorPane;
    public FXMLDocumentController(){
        
    }
    
    int currentX, currentY, currentZ;

    //FXML Link to quit button    
    @FXML
    private void byebye(ActionEvent event) {
        
        ((Node)(event.getSource())).getScene().getWindow().hide();
        
    }
    
  
    @FXML
    private void clearHandle(ActionEvent event) {
        
        
       
        sliderX.setValue(0);
        sliderY.setValue(0);
        sliderZ.setValue(0);
        sliderT.setValue(0);
        sliderSR.setValue(0);
        sliderSimTime.setValue(0);
        
    }
    

    //FMXL Link to run button
    @FXML
    private void runSim(ActionEvent event) {
        
        heartViewer.setImage(heartbeating);
        
        nx = Integer.parseInt(xLabel.getText());
        ny = Integer.parseInt(yLabel.getText());
        nz = Integer.parseInt(zLabel.getText());
        
        thresh = Integer.parseInt(threshLabel.getText());
        
        
        try {
                    
            PrintWriter out = new PrintWriter("/home/varderes/Desktop/CellSimulation/params.txt");
            out.println(nx);
            out.println(ny);
            out.println(nz);
            out.println(thresh);
            out.close();
                    
            } catch (IOException e) {
             
        }
        
        
        copyWorker = createWorker();
        
        progressBar.progressProperty().unbind();
        progressBar.progressProperty().bind(copyWorker.progressProperty());
        
        new Thread(copyWorker).start();
        
        System.out.println(nx + " " + ny + " " + nz + " " + thresh);
        
    }
    
    
    
    @Override
    public void initialize(URL url, ResourceBundle rb) {
        // TODO
        
        selection.setItems(selections);
        
        
        
        
        
        sliderX.valueProperty().addListener(new ChangeListener<Number>() {
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                xLabel.setText(String.valueOf((newValue.intValue())));
                
                currentX = newValue.intValue();
                
                if(currentX > 3 && currentY > 3 && currentZ > 3) {
                    run.setDisable(false);
                    warning1.setFill(Color.GREEN);
                    warning1.setText("Lengths set!");
                }
                else {
                    run.setDisable(true);
                    warning1.setFill(Color.RED);
                    warning1.setText("Minimum lengths must be 3!");
                }
                     
              
            }
        });
        
        
        sliderY.valueProperty().addListener(new ChangeListener<Number>() {
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                yLabel.setText(String.valueOf((newValue.intValue())));
                
                currentY = newValue.intValue();
                
                if(currentX > 3 && currentY > 3 && currentZ > 3) {
                    run.setDisable(false);
                    warning1.setFill(Color.GREEN);
                    warning1.setText("Lengths set!");
                }
                else {
                    run.setDisable(true);
                    warning1.setFill(Color.RED);
                    warning1.setText("Minimum lengths must be 3!");
                }
                
            }
        });
        
        
        sliderZ.valueProperty().addListener(new ChangeListener<Number>() {
            
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                zLabel.setText(String.valueOf((newValue.intValue())));
                
                
                currentZ = newValue.intValue();
                
                if(currentX > 3 && currentY > 3 && currentZ > 3) {
                    run.setDisable(false);
                    warning1.setFill(Color.GREEN);
                    warning1.setText("Lengths set!");
                }
                else {
                    run.setDisable(true);
                    warning1.setFill(Color.RED);
                    warning1.setText("Minimum lengths must be 3!");
                }
                
            }
            
        });
        
        
        sliderT.valueProperty().addListener(new ChangeListener<Number>() {
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                threshLabel.setText(String.valueOf((newValue.intValue())));
            }
        });
        
        
        sliderSimTime.valueProperty().addListener(new ChangeListener<Number>() {
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                simtimeLabel.setText(String.valueOf((newValue.intValue()))+ "ms");
            }
        });
        
        
        sliderSR.valueProperty().addListener(new ChangeListener<Number>() {
        @Override
            public void changed(ObservableValue<? extends Number> observable,
                Number oldValue, Number newValue) {

                srloadLabel.setText(String.valueOf((newValue.intValue())));
            }
        });
        
        
        
        selection.getSelectionModel().selectedIndexProperty().addListener(
            new ChangeListener<Number>() {
                public void changed(ObservableValue ov, Number value, Number new_value) {
                    
                    if(new_value != null) {
                        optionWarning.setFill(Color.GREEN);
                        optionWarning.setText("Option Selected!");
                    }
                    
                }
            });
        
        
        
    }
    
    

    public Task createWorker() {
    return new Task() {
        @Override

        protected Object call() throws Exception {

            // Running C++ simulation
            String s = null;
            String s2 = null;
            String s3 = null;
            String[] env = {"/home/varderes/Desktop/CellSimulation"};
            String cmd = "./runcpp.sh";
            String cmd2 = "./run.sh";


            try {

                Process p = Runtime.getRuntime().exec("/home/varderes/Desktop/CellSimulation/runcpp.sh");
                
                
                BufferedReader stdInput = new BufferedReader(new
                     InputStreamReader(p.getInputStream()));
                BufferedReader stdError = new BufferedReader(new
                     InputStreamReader(p.getErrorStream()));
                
                
                System.out.println("Here is the standard output of the command:\n");

                while ((s = stdInput.readLine()) != null) {
                    if(s != null && s.length() > 0) {
                        try {

                            double value = Double.parseDouble(s)/100;

                            status.setText("Running Simulation");
                            updateProgress(value, 1);
                            progressText.setText(String.valueOf((int)(value*100)) + "%");

                            if(value == 1.0){
                                status.setText("Simulation Complete");
                                status.setFill(Color.GREEN);
                                
                            }

                            if(value == 2.0){
                                status.setText("Importing Data to Mathematica");
                                status.setFill(Color.DARKVIOLET);
                            }

                            if(value == 3.0){
                                status.setText("Creating Plots");
                            }

                            if(value == 4.0){
                                status.setText("Saving Plots");
                            }

                            if(value == 5.0){
                                status.setText("Opening Mathematica!");
                                status.setFill(Color.GREEN);
                            }

                        } catch (Exception e) {

                        }
                    }
                }




                // read any errors from the attempted command
                System.out.println("Here is the standard error of the command:\n");
                while ((s = stdError.readLine()) != null) {
                    System.out.println(s);
                }     
            } catch (IOException e) {
                System.out.println("exception happened - here's what I know: ");
                e.printStackTrace();
                //System.exit(-1);
            }                


            
            try {
                

                Process p2 = Runtime.getRuntime().exec("/home/varderes/Desktop/CellSimulation/run.sh");
                
                BufferedReader stdInput = new BufferedReader(new
                     InputStreamReader(p2.getInputStream()));
                BufferedReader stdError = new BufferedReader(new
                     InputStreamReader(p2.getErrorStream()));
                
                
                System.out.println("Here is the standard output of the command:\n");
                while ((s2 = stdInput.readLine()) != null) {
                    if(s2 != null && s2.length() > 0) {


                        System.out.println(s2);


                        try {

                            double value2 = Double.parseDouble(s2);

                            if(value2 == 2.0){
                                status.setText("Importing Data to Mathematica...");
                                status.setFill(Color.DARKVIOLET);
                            }

                            if(value2 == 3.0){
                                status.setText("Creating Plots...");
                            }

                            if(value2 == 4.0){
                                status.setText("Saving Plots...");
                            }

                            if(value2 == 5.0){
                                status.setText("Opening Mathematica!");
                                status.setFill(Color.GREEN);
                            }

                        } catch (Exception e) {

                        }

                    }

                }
                // read any errors from the attempted command
                System.out.println("Here is the standard error of the command:\n");
                while ((s2 = stdError.readLine()) != null) {
                    System.out.println(s2);
                }     
            } catch (IOException e) {
                System.out.println("exception happened - here's what I know: ");
                e.printStackTrace();
                //System.exit(-1);
            }    


            try {

                Process p3 = Runtime.getRuntime().exec("/home/varderes/Desktop/CellSimulation/mathematica.sh");
                BufferedReader stdInput = new BufferedReader(new
                     InputStreamReader(p3.getInputStream()));
                BufferedReader stdError = new BufferedReader(new
                     InputStreamReader(p3.getErrorStream()));
                // read the output from the command
                System.out.println("Here is the standard output of the command:\n");
                while ((s3 = stdInput.readLine()) != null) {
                    System.out.println(s3);

                }
                // read any errors from the attempted command
                System.out.println("Here is the standard error of the command:\n");
                while ((s3 = stdError.readLine()) != null) {
                    System.out.println(s3);
                }     
            } catch (IOException e) {
                System.out.println("exception happened - here's what I know: ");
                e.printStackTrace();
                System.exit(-1);
            }
            return true;
        }
    };                    
                
    }               
	
    
    
    
    
}



