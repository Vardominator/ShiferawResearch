/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cellsimfxml;

import javafx.application.Application;
import static javafx.application.Application.launch;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

/**
 *
 * @author varderes
 */
public class CellSimFXML extends Application {
    
        
    @Override
    public void start(Stage mainStage) throws Exception {
     
        
        // Link FXML file
        Parent root = FXMLLoader.load(getClass().getResource("design.fxml"));
        
        Scene scene = new Scene(root);
        
        // Link CSS file
        scene.getStylesheets().add(CellSimFXML.class.getResource("CascadeStyleSheet.css").toExternalForm());
                
        mainStage.setResizable(false);
        mainStage.setTitle("Heart Cell Simulator");
        mainStage.setScene(scene);
        mainStage.show();

    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        launch(args);
        
    }
    
}
