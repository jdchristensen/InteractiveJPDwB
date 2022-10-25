// Interactive demo of the persistence pipeline in 2D,
// building on javaplexDemo.pde from https://github.com/appliedtopology/javaplex.

// To do:
// - drawing the simplicies is the bottle-neck, not doing the homology calculations.
//   (Can see this using the 'h' keystroke, since then the homology is not computed.)
//   How to optimize?
// - If f is increased, keep existing image and draw only the additional simplices.
//   Should be much faster.  No such trick works for decreasing f, though.

import edu.stanford.math.plex4.api.*;
import edu.stanford.math.plex4.examples.*;
import edu.stanford.math.plex4.streams.impl.VietorisRipsStream;
import edu.stanford.math.plex4.homology.chain_basis.Simplex;
import edu.stanford.math.plex4.homology.filtration.FiltrationConverter;
import edu.stanford.math.plex4.homology.interfaces.AbstractPersistenceAlgorithm;
import edu.stanford.math.plex4.homology.barcodes.*;
import java.util.Map.Entry;
import java.util.List;

double[][] pts;
double[][] pts1, pts2, pts3, pts4, pts5, pts6, pts7, pts8;

float sizeX=7,sizeY=7;
double eps;
double f;
double maxeps; // maximum epsilon value passed to javaplex
VietorisRipsStream<double[]> vrs;
FiltrationConverter fc;
AbstractPersistenceAlgorithm<Simplex> algo;
BarcodeCollection<Double> ints=null;
float[][] intervals;
int num_pts = 0;
PFont ft;
float int_max = 0;
float max = 0;
int barcode_width;
int barcode_height;
float barcode_left;
float barcode_top;
int font_size = 20;
boolean show_barcode = true;
boolean redraw_needed = true;
boolean recalc_needed = true;
int cur_width = 0, cur_height = 0;

void settings() {
  // It works best to start at the right size initially.
  // But fullScreen() doesn't get it quite right.
  //fullScreen();
  size(2560, 1440);
  //size(1920, 1050);
}

void setup() {
  // Most things adapt to the window resizing, but the loaded point cloud
  // data is scaled to the starting window size...  
  surface.setResizable(true);

  ft = createFont("Courier", 26, true);
  textFont(ft, font_size);
  resetPoints();
  pts1 = load_pts(1);
  pts2 = load_pts(2);
  pts3 = load_pts(3);
  pts4 = load_pts(4);  
  pts5 = load_pts(0);
  pts6 = load_pts(0);
  pts7 = load_pts(0);
  pts8 = load_pts(0);
}

void draw() {
  if(recalc_needed)
  {
    setupVRS();    
    recalc_needed = false;
  }

  // Detect resizing:
  if((cur_width != width) || (cur_height != height))
  {
    redraw_needed = true;
    cur_width = width;
    cur_height = height;
  }

  if(!redraw_needed)
  {
    delay(50); // milliseconds
    return;
  }

  barcode_width = (int)(0.42*width);
  barcode_height = (int)(0.75*height);
  // Slightly off center to the right, to make room for "dim i" labels:
  barcode_left = 0.75*width - barcode_width/2 + font_size*2;
  // Slightly up from center to make room for the x-axis labels:
  barcode_top = 0.4*height - barcode_height/2 - font_size;

  background(255);

  line(width/2.0, 0, width/2.0, 0.8*height);
  line(0, 0.8*height, width, 0.8*height);

  draw_instructions(0.05*width, 0.8*height, width, height);
  draw_barcode();

  strokeWeight(2);

  for(Simplex s : vrs) {
    double fv = fc.getFiltrationValue(vrs.getFiltrationIndex(s));

    if(fv > f)
      continue;

    int[] ix;
    ix = s.getVertices();

    switch(s.getDimension()) {
      case 0:
        fill(0);
        ellipse((float)pts[ix[0]][0],(float)pts[ix[0]][1],sizeX,sizeY);
        break;
      case 1:
        fill(0);
        line((float)pts[ix[0]][0],(float)pts[ix[0]][1],
            (float)pts[ix[1]][0],(float)pts[ix[1]][1]);
        break;
      case 2:
        fill(0,0,255,30);
        triangle((float)pts[ix[0]][0],(float)pts[ix[0]][1],
            (float)pts[ix[1]][0],(float)pts[ix[1]][1],
            (float)pts[ix[2]][0],(float)pts[ix[2]][1]);
        break;
      default:
        continue;
    }
  }
  strokeWeight(1);
  redraw_needed = false;
}

//*****************************************
// Compute a new VietorisRipsStream
//*****************************************

void setupVRS() {
  // 2 is maxDimension; maxeps is maxFiltrationValue, 1000 is numDivisions.
  // Putting numDivisions lower doesn't speed things up noticably.
  vrs = Plex4.createVietorisRipsStream(pts,2,maxeps,1000);
  fc = vrs.getConverter();

  redraw_needed = true;
  
  // We'll compute the barcodes even when not displaying them,
  // since this is used to know the maximum epsilon value.
  //if(!show_barcode) return;

  // compute intervals
  // "2" means to compute H_0 and H_1.  I'm guessing this defaults
  // to working over the rationals?
  algo = Plex4.getDefaultSimplicialAlgorithm(2);
  ints = algo.computeIntervals(vrs);
  //println(ints);
  String s = ints.toString();
  
  // convert intervals into tidy array
  intervals = ints_to_intervals(s);
}

//*****************************************
// Reset the points buffer
//*****************************************

void resetPoints() {
      pts=new double[0][2];
      f = 0;
      eps = 10;
      maxeps = 400; // For smaller data sets, ~1000 works ok.
}

//*****************************************
// Display instructions at bottom
//*****************************************

void draw_instructions(float xa, float ya, float xb, float yb) {
  int h = font_size-1;
  fill(0);
  //text("INSTRUCTIONS", xa+10, yb-10*h);
  text("click in fig -- adds a point (shift-click to remove)", xa +10, yb-9*h);
  text("click in VR  -- set Vietoris-Rips parameter", xa +10, yb-8*h);
  text("1-4          -- loads example data sets", xa+10, yb-7*h);
  text("T, Y, U, I   -- saves current data set to 5, 6, 7, or 8", xa+10, yb-6*h);
  text("5-8          -- loads saved data set", xa+10, yb-5*h);  
  text("C            -- clear points", xa + 10, yb-2*h);
  text("RIGHT/LEFT   -- step Vietoris-Rips parameter forward/back", xa+10, yb-4*h);
  text("0            -- set Vietoris-Rips parameter to 0", xa+10, yb-3*h);
  text("Q            -- quit", xa+10, yb-h);
  //text("InteractiveJPDwB Instructions.", xa+600, yb-9*h);
  text("Luke Wolcott, 2016, Dan Christensen, 2018.", xa+600, yb-h);
}

//*****************************************
// On mousepress, if within data box then add a point.
// If in barcode, set f based on x position.
//*****************************************

void mousePressed() {
  if (keyPressed && keyCode == SHIFT){ 
    for (int i = 0; i < pts.length; i++){
      if (sq((float) pts[i][0]-mouseX)+sq((float) pts[i][1]-mouseY) < 25){  
        pts = remove_row(i);
        recalc_needed = true;
        break;
      }
    }
  }
  else{ 
    if ((mouseX < width/2) && (mouseY < 0.8*height)) {
      double[] pt = new double[2];

      pt[0] = mouseX;
      pt[1] = mouseY;
      
      println(pt[0]+","+pt[1]);
      pts = (double[][]) append(pts,pt);
      recalc_needed = true;
    }
    else if (show_barcode && (mouseX >= barcode_left) && (mouseX <= barcode_left+barcode_width)
             && (mouseY >= barcode_top) && (mouseY <= barcode_top+barcode_height+40)) {
      f = ((mouseX - barcode_left)/barcode_width)*max;
      redraw_needed = true;
    }
  }
}

//*****************************************
// On keypress:
//
// Q      -- quit
// C      -- clear points
// LEFT   -- step Vietoris-Rips complex back
// RIGHT  -- step Vietoris-Rips complex forward
// 1-4    -- loads pre-stored data sets
// 5-8    -- loads data sets saved by the user
// T,Y,U,I-- saves data set in 5-8, respectively
//*****************************************

void keyPressed() {
  switch(key) {
    case 'q':
    case 'Q':
      exit();
      break;
      
    case 'c':
    case 'C':
      resetPoints();
      recalc_needed = true;
      break;
      
    case '1':                              
      pts = pts1;
      recalc_needed = true;
      break;    
  
    case '2':                              
      pts = pts2;
      recalc_needed = true;
      break;
     
    case '3':                              
      pts = pts3;
      recalc_needed = true;
      break;
     
    case '4':                             
      pts = pts4;
      recalc_needed = true;
      break;
      
    case '5':
      pts = pts5;
      recalc_needed = true;
      break;
      
    case '6':
      pts = pts6;
      recalc_needed = true;
      break;
      
    case '7':
      pts = pts7;
      recalc_needed = true;
      break;
      
    case '8':
      pts = pts8;
      recalc_needed = true;
      break;
      
    case 't':
    case 'T':
      pts5 = pts;
      break;
    
    case 'y':
    case 'Y':
      pts6 = pts;
      break;
    
    case 'u':
    case 'U':
      pts7 = pts;
      break;
    
    case 'i':
    case 'I':
      pts8 = pts;
      break;    
         
    case '0':
      f = 0;
      redraw_needed = true;
      break;
      
    case 'h':
      show_barcode = !show_barcode;
      redraw_needed = true;
      break;

    case CODED:
      switch(keyCode) {
        case RIGHT:
          if(f < 200)
            f += eps;
          else
            f += 2*eps;
          // Using maxeps instead of max here means that f can go beyond
          // the barcode plot.  Sometimes that's useful.
          if (f>maxeps)
            f=maxeps;
          redraw_needed = true;
          println(f+": "+eps);
          break;
        case LEFT:
          if(f < 200)
            f -= eps;
          else
            f -= 2*eps;
          if(f<0)
            f=0;
          redraw_needed = true;
          println(f+": "+eps);
          break;
    }    
  }
}

//*****************************************
// Load pre-saved data into tables so it is ready to be read with 1-4 and 5-8.
//*****************************************

double[][] load_pts(int n) {
  double[][] ptsn;
  Table tablen;
  if (n==0){
    ptsn = new double[0][2];
    return ptsn;
  }
  else if (n==1)
    tablen = loadTable("seed1_data.csv", "header");
  else if (n==2)
    tablen = loadTable("seed2_cross.csv", "header");
  else if (n==3)
    tablen = loadTable("seed3_cross_hole_005.csv", "header");
  else 
    tablen = loadTable("seed4_circle.csv", "header");              
  ptsn = new double[tablen.getRowCount()][2];
  int i=0;
  for (TableRow row : tablen.rows()){
    ptsn[i][0] = row.getDouble("X_value")*width/2*0.9+width/2*0.05;
    ptsn[i][1] = row.getDouble("Y_value")*0.8*height*0.9+0.8*height*0.05;
    i++;
  }
  return ptsn;
}

//*****************************************
// Convert the output of a Javaplex persistence interval calculation into a tidy array.
//*****************************************

float[][] ints_to_intervals(String s){    
  String[] lines = splitTokens(s, " [,)\n\r");  

//  print how many lines, and then print each line
//  println("there are " + lines.length + " lines");
//  for (int i=0; i<lines.length; i++){
//    println(lines[i]);
//  }

  // Count how many different dimensions and points.
  int num_dims = 0;
  for (int i=0; i<lines.length; i++){
    if (lines[i].equals("Dimension:")){
      num_dims = num_dims + 1;
    }
  }      

  // print how many dimensions there are, and how many intervals there are
  println("there are " + num_dims + " different dimensions");
  num_pts = lines.length/2 - num_dims;
  println("there are " + num_pts + " different intervals");            

  // Build array of dimension, start, end.
  intervals = new float[num_pts][3];
  int dim = 0;
  int pt_number=-1;
  for (int k = 0; k<lines.length/2; k++){
    if (lines[2*k].equals("Dimension:")){
      dim = int(lines[2*k+1]);
      //println("dim:" + dim);
    }
    else{
      pt_number = pt_number + 1;
      //println(pt_number);
      intervals[pt_number][0] = dim;
      intervals[pt_number][1] = float(lines[2*k]);
      intervals[pt_number][2] = float(lines[2*k+1]);
    }  
  }  

  // print array of intervals
  for( int j=0; j<num_pts; j++){
    println(intervals[j][0], intervals[j][1], intervals[j][2]);
  } 
  
  // Figure out horizontal scale.
  int_max = 0;
  for (int i=0; i<intervals.length; i=i+1){
    if (intervals[i][2] > int_max){
      int_max = intervals[i][2];
    }
    // The above ignore's NaN's, so also need to handle when the start of
    // an infinite interval is larger than the ends of all of the finite
    // intervals:
    if (intervals[i][1] > int_max){
      int_max = intervals[i][1];
    }
  }
  println("Max interval value is " + int_max + ".");
  max = int_max+100;           // Adjust so max isn't cut off.
  if(max < 400) max = 400;     // Ensure at least some range given.

  return intervals;
}


//*****************************************
// Draw the barcode corresponding to a tidy array of intervals.
//*****************************************

void array_to_barcode(float[][] intervals){
  int nrow = intervals.length;
       
  // Look through table and figure out where the dimension changes.  
  int[] spots = {0};
  for (int i=1; i<nrow; i = i + 1){
    if (intervals[i][0] > intervals[i-1][0] ){
      spots = (int[]) append(spots,i);
    }
  }
    
  // Figure out what those dimensions actually are.  
  int[] dims = {int(intervals[0][0])};
  for (int i=1; i<spots.length; i=i+1){
    dims = (int[]) append(dims, int(intervals[spots[i]][0]));
  }
  spots = (int[]) append(spots,nrow);
  
  // Start drawing lines.
  // These are set in draw() so they can be used for mouse presses too:
  float a = barcode_left;
  float b = barcode_top;
  float spaces = nrow - dims.length + 2 + 2 +4*(dims.length-1);
  float incr = barcode_height/spaces;
  float y=b;      
  textSize(font_size);
  fill(0);
  for (int j=0; j<dims.length; j=j+1){
    y = y + 2*incr;
    text("dim " + dims[j], a-font_size*4, y);
    if(dims[j] == 1)
      stroke(255,0,0);
    for (int k=spots[j]; k<spots[j+1]; k=k+1){    
      float start = intervals[k][1];
      float finish = intervals[k][2];
      // Convert infinity to max length.
      if(Float.isNaN(finish)) finish = max;
      strokeWeight(2.5);
      line(a + barcode_width*(start/max), y, a + barcode_width*(finish/max), y);
      y = y + incr;
    }
    if (j < (dims.length-1)){
      y = y + 1*incr;
      stroke(0,0,204);
      strokeWeight(1.5);
      line(a, y, a+barcode_width, y);        // Draws a full line to separate dimensions
      stroke(0);
      strokeWeight(1);
    }
  }
  text(int(max), a+barcode_width-30, b+barcode_height+font_size+10);
  text(0, a-2, b+barcode_height+font_size+10);
  stroke(0,0,204);
  strokeWeight(2);      
  line(a, b+barcode_height, a, b+barcode_height+15);
  line(a+barcode_width, b+barcode_height, a+barcode_width, b+barcode_height+15);
  line(a + barcode_width*((float) f/max), b + barcode_height,
       a + barcode_width*((float) f/max), b + barcode_height+15);
  stroke(0);
  strokeWeight(1);
}
  
void draw_barcode(){
      //initialize barcode region   
      fill(255);
      rect(width/2, 0, width, 0.8*height);
      
      if(!show_barcode) return;
      
      stroke(0,0,204);
      strokeWeight(1.5);
      rect(0.75*width - barcode_width/2 + font_size*2, 0.4*height - barcode_height/2 - font_size, barcode_width, barcode_height);
      stroke(0);
      strokeWeight(1);
      
      // quit if there are no points in the array
      if (intervals.length == 0) {
        return;
      }
      
      // draw barcode from tidy array
      array_to_barcode(intervals);
}
  
//*****************************************
// Remove a row from pts.
//*****************************************
  
double[][] remove_row(int r){
  double[][] temp = new double[0][2];
  for (int i=0; i<pts.length; i++){
    if (i != r)
      temp = (double[][]) append(temp, pts[i]);
  }
  return temp;
}
