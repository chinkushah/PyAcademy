import math
print("Welcome to Py Academy")

for ch in range(1,100):

  what_subject = input("What subject may I help you with today?(Math/Physics/Chemistry) ")

    
  if what_subject == "math" or what_subject == "Math" or what_subject == "mathematics" or what_subject == "Mathematics":

    what_chap = input("What chapter may I help you with?(Progressions/Straight Lines(sl)/ Calculus)")

    if what_chap == "progressions" or what_chap == "Progressions":

      print("The topics involved along with their formulae are:")
      print('''For any Arithmetic Progression;
      a = first term of an AP, d = common difference of an AP
      nth term of an AP: a + (n-1)d
      Sum of n terms of an AP: n/2[2a + (n-1)d]

        Arithmetic Mean of two numbers 'a' and 'b';
        AM = [a+b]/2
        d = [b-a]/[n+1]

        For any Geometric Progression;
        a = first term of a GP, r = common ratio of GP
        nth term of a GP: ar^(n-1)
        Sum of n terms of GP: [a(r^n - 1)]/r-1

        Geometric Mean of two numbers 'a' and 'b';
        GM = (ab)^[1/2]
        r = (b/a)^[1/n+1]''')

      more_help=input("Do you need further assistance?(Yes/No) ")

      if more_help == "yes" or more_help == "Yes":
        ProgOrMean = (input("Do you want to find AM/GM or nth term/sum or insertion of AM/GM(ins)? "))

        if ProgOrMean == "nth term/sum" or ProgOrMean == "nth term/Sum":

          first_term = input("Enter First Term of the Progression: ")
          first_term = float(first_term)
          is_ap_or_gp = input("Is the Progression AP or GP?")
          is_ap_or_gp = str(is_ap_or_gp)

          if is_ap_or_gp == "AP" or is_ap_or_gp == "ap":
            common_difference = input("Enter Common Difference:")
            common_difference = float(common_difference)
            term = input("Enter the Term:")
            term = int(term)
            find_nth_term_or_sum = input("Do You Want to Find nth term or sum? ")
            find_nth_term_or_sum = str(find_nth_term_or_sum)

            if find_nth_term_or_sum == "nth term" or find_nth_term_or_sum == "nth Term":
              nth_term = first_term + ((term - 1) * common_difference)
              print("the nth Term is", nth_term)

            elif find_nth_term_or_sum == "sum" or find_nth_term_or_sum == "Sum":
              Sum = (term/2)*((2*first_term) + ((term-1)*common_difference))     
              print("The Sum is", Sum)
          else:
            common_ratio = input("Enter Common Ratio of GP:" )
            common_ratio = float(common_ratio)
            term = input("Enter nth Term of GP:")
            term = int(term)
            find_nth_term_or_sum = input("Do You Want to Find nth term or sum?")

            if find_nth_term_or_sum == "nth term" or find_nth_term_or_sum == "nth Term":
              nth_term = round(((first_term)*((common_ratio)**(term-1)),2)) 
              print("The nth Term is", nth_term)

            elif find_nth_term_or_sum == "sum" or find_nth_term_or_sum == "Sum":
              Sum = ((first_term*(1-common_ratio**term))/(1-common_ratio))   
              print("The Sum is", Sum)

        elif ProgOrMean == "AM/GM" or ProgOrMean == "am/gm":
          AM_GM = input("Do you want to find AM or GM?")

          if AM_GM == "AM" or AM_GM == "am":
            term_one = int(input("Enter one term:"))
            term_two = int(input("Enter second term:"))
            AM = (term_one + term_two)/2
            print("The AM is",AM)

          else:
            term_one = int(input("Enter one term:"))
            term_two = int(input("Enter second term:"))
            GM = (term_one*term_two)**(1/2)
            print("The GM is",GM)

        else:

          AMorGM = input("Insertion of AMs or GMs?")

          if AMorGM == "AM" or AMorGM == "AMs":
            a = int(input("Enter first term: "))
            b = int(input("Enter last term: "))
            n = int(input("Enter the number of terms you want to enter: "))
            d = (b-a)/(n+1)
            series = 0
            print("The AP thus formed is")

            for ch in range(0,n+2):
              Series = a + (d*ch)
              print(Series)

          else:
            
            a = int(input("Enter first term: "))
            b = int(input("Enter last term: "))
            n = int(input("Enter the number of terms you want to insert: "))
            r = (b/a)**(1/(n+1))
            series = 1
            print("The GP thus formed is")

            for ch in range(0,n+2):
              Series = a*(r**ch)
              print(Series)
              
    
    
    elif what_chap == 'straight lines' or what_chap == 'sl':

      print('''The topics involved along with their formulae are:
      General equation of a line is ax + by + c = 0.

      If equation of a line is of the form y = mx+c, then m is the slope of the line.


        Slope of a line given two points (a,b) and (c,d);
        (d-b)/(c-a) = (y-b)/(x-a).

        Angle(A) between two lines with slopes m and M ;
        tanA = (M-m)/(1+mM).''')

      more_help = input("Do you need further assistance?")

      if more_help == "yes" or more_help == "Yes":
                
        dist = input("Do you want to find the distance of a point from a line?")

        if dist == "yes" or dist == "Yes":

          y_coordinate = float(input("Enter y-coordinate of the point:"))
          x_coordinate = float(input("Enter x-coordinate of the point:"))
          coeff_y = float(input("Enter coefficient of y from the equation of the line:"))
          coeff_x = float(input("Enter coefficient of x from the equation of the line:"))
          constant = float(input("Enter constant term from the equation of the line:"))
          distance = round((y_coordinate*coeff_y + x_coordinate*coeff_x + constant)/((coeff_x**2) + (coeff_y**2))**(1/2),2)
          print("The perpendicular distance of the point from the line is",distance)

        else:
          coordinates_given = input("Are the coordinates of line given?")

          if coordinates_given == "yes" or coordinates_given == "Yes":
            y1 = float(input("Enter first y-coordinate:"))
            y2 = float(input("Enter second y-coordinate:"))
            x1 = float(input("Enter first x-coordinate:"))          
            x2 = float(input("Enter second x-coordinate:"))
            slope = ((y2-y1)/(x2-x1))
            print("The slope of the line is",slope)
            y_diff = y2-y1
            x_diff = x1-x2
            constant = (x1*(y1-y2) + y1*(x2-x1))
            angle = round(math.degrees(math.atan(slope)),2)
            print("The angle made by the line with the x-axis is",angle,"degrees")
            print("The equation of the line is",y_diff,"x +",x_diff,"y" "+",constant,"= 0")
            
            from matplotlib import pyplot as plt
            plt.plot([x1,x2],[y1,y2])
            plt.show()
                    
          else:
              
            slope = float(input("Enter slope of the line:"))
            y_int = float(input("Enter y-intercept of the line:"))
            print("The equation of the line is y =", slope,"x +", y_int)
            from matplotlib import pyplot as plt
            plt.plot([0,(-y_int/slope)],[y_int,0])
            plt.show()
            
    elif what_chap == 'c' or what_chap == 'Calculus':
      from sympy import *
      import matplotlib.pyplot as plt
      x = Symbol('x')
      y = Symbol('y')
      calc = input("Do you want to differentiate or integrate a function? (diff/int)")

      if calc == 'diff':
          f = input("Enter function to be differentiated :")
          print(diff(f,x))
              
      else:
          f = input("Enter function to be integrated :")
          print(integrate(f,x))
      continue
      


  elif what_subject == "physics" or what_subject == "Physics":
    what_chap = input("What chapter do you need help with, Projectile Motion(pm) or Circular Motion(cm)? ")

    if what_chap == "projectile motion" or what_chap == "Projectile Motion" or what_chap == "Projectile motion" or what_chap == "pm":

      x = float(input("Enter Initial Velocity(m/s):"))
      t = float(input("Enter Angle of Projection(degrees):"))
      y = math.radians(t)

      time_of_flight = ((x*(math.sin(y)))/5)
      print("Time of Flight is",time_of_flight,"seconds")

      horizontal_range = (((x**2)*(math.cos(y))*(math.sin(y)))/5)
      print("Horizontal Range of the Projectile is",horizontal_range,"meters")

      maximum_height = (((x**2)*((math.sin(y)**2)))/20)
      print("Maximum Height of the Projectile is",maximum_height,"meters")

      coeff_x = (5/(x*math.cos(y))**2)
      eqn = ('y =',math.tan(y),'x -',coeff_x,'x^2')
      print("The equation of the projectile is")
      print('y =',math.tan(y),'x -',coeff_x,'x^2')

    elif what_chap == "Circular Motion" or what_chap == "circular motion" or what_chap == "cm":
      find = input("What do you want to find, Angular Velocity(av), Angular Acceleration(aa)? ")

      if find == "angular velocity" or find == "Angular Velocity" or find == "av":
        accn_giv = input("Is the angular acceleration given?")

        if accn_giv == "Yes" or accn_giv == "yes":
          ang_accn = float(input("Enter the angular acceleration(in rad/s^2):"))
          ang_disp = float(input("Enter the angular displacement(in rad):"))
          ang_vel = (2*ang_accn*ang_disp)**(1/2)
          print("The angular velocity is",ang_vel,"rad/s")

        else:
          cent_accn = input("Is the centripetal acceleration given?")

          if cent_accn == "yes" or cent_accn == "Yes":
            cent_accn == float(input("Enter the centripetal acceleration(in m/s^2):"))
            radius = float(input("Enter the radius of circular motion(in m):"))
            vel = (cent_accn*radius)**(1/2)
            ang_vel = (vel/radius)
            print("The angular velocity is",ang_vel,"rad/s")

          else:
            lin_accn = float(input("Enter the linear acceleration(in m/s^2):"))
            radius = float(input("Enter the radius of circular motion(in m):"))
            ang_disp = float(input("Enter the angular displacement(in rad):"))
            ang_accn = lin_accn/radius
            ang_vel = (2*ang_accn*ang_disp)**(1/2)
            print("The angular velocity is",ang_vel,"rad/s")

      elif find == "angular acceleration" or find == "Angular Acceleration" or find == "aa":
        ang_vel = input("Is the angular velocity given?")

        if ang_vel == "Yes" or ang_vel == "yes":
          ang_vel = float(input("Enter the angular velocity(in rad/s):"))
          ang_disp = float(input("Enter the angular displacement(in rad):"))
          ang_accn = (ang_vel/(2*ang_disp))
            
        else:
          cent_accn = input("Is the centripetal acceleration given?")

          if cent_accn == "Yes" or cent_accn == "yes":
            cent_accn = float(input("Enter the centripetal acceleration(in m/s):"))
            accn = float(input("Enter net acceleration(in m/s^2):"))
            ang_accn = ((accn**2)-(cent_accn**2))**(1/2)
            print("The angular acceleration is",ang_accn)


  elif what_subject == "Chemistry" or what_subject == "chemistry":
    import pandas as pd   
    df = pd.read_csv("") #ENTER 'PERIODIC TABLE OF ELEMENTS.csv' FILE LOCATION IN THE QUOTES TO THE LEFT
    df = pd.DataFrame(df)
    df["MassNumber"] = df["NumberofNeutrons"] + df["NumberofProtons"]
    df_atno = df.set_index("AtomicNumber")
    df_atmass = df.set_index("AtomicMass")
    df_massno = df.set_index("MassNumber")
    

    chem_chp = input("Which chapter do you need help with? Mole Concepts (mc) or Atomic Structure (as): ")

    if chem_chp == "Mole concepts" or chem_chp == "mole concepts" or chem_chp == "Mole Concepts" or chem_chp == "mc":
        
        print('''Here are some helpful formulae of this chapter:
        No. of moles: n = (Given mass)/(Molar mass)
                = (No. of particles)/Na                   _Where Na is Avogadros's Nomber.
                = (Vol. of gas at STP)/(22.4l)

        Average Atomic Mass of Elements having different isotopes: Mavg= (M1a1 + M2a2 + M3a3..Mnan)/(a1 + a2 + a3...+an)   _where a is percentage of abundance of isotope.
        Mass percent of an element = (mass of that element in compound*100)/(mass of compound)

        Vapour density: (Molar Mass of gas)/2

        Molarity: M = (moles of solute)/(volume of solution)
        Molality: m = (moles of solute)/(mass of solvent)Z
        Mole fraction: Of solute A (XA) = nA/(nA + nB) , Of solvent B (XB) = nB/(nA + nB)''')

        mole_concept_notes = input("Do you require any further assisstance? ")
        if mole_concept_notes == "Yes" or mole_concept_notes == "yes":
          
            help_mole_concept = input("What do you need help with? Mass Percent (mp) , Molarity , Molality , Empirical Formula (ef)  ")

            if help_mole_concept == "Mass Percent" or help_mole_concept == "mass percent" or help_mole_concept == "mp":            
                totalMass = 0
                elements = int(input("How many elements are present in the compound?"))

                for mass in range(1,elements + 1):
                    Atmass = input("Enter the element:")
                    atomicMass = float(df_atmass[df_atmass["Symbol"] == Atmass].index.values)
                    NumMolecule = int(input("Enter the number of atoms of the particular element: "))
                    mass =  atomicMass * NumMolecule
                    totalMass += mass
                print("The mass of this compound is",totalMass)

                Element = input("Which element's mass percent would you like to find? ")
                moles = float(input("Give number of atoms of element: "))
                Mass = float(df_atmass[df_atmass["Symbol"] == Element].index.values*moles)
                print("Mass of element is atomic mass*moles = ", Mass)
                print("Mass Percent of the element is: ", Mass*100/totalMass)

            elif help_mole_concept == "Molarity" or help_mole_concept == "molarity":
                moles = float(input("Give moles of element: "))
                vol = float(input("Give volume of solution: "))
                print("Molarity =", moles/vol )

            elif help_mole_concept == "Molality" or help_mole_concept == "molality":
                moles = float(input("Give moles of element: "))
                mass = float(input("Give mass of solvent in kg: "))
                print("Molality= ", moles/mass)

            elif help_mole_concept == "Empirical Formula" or help_mole_concept == "empirical formula" or help_mole_concept == "ef":
                totalMass = 0
                elements = int(input("How many elements are present in the compound?"))

                if elements == 3:

                    ele1 = input("Enter the element: ")
                    per1 = float(input("Percentage of this element: "))
                    ele2 = input("Enter the element: ")
                    per2 = float(input("Percentage of this element: "))
                    ele3 = input("Enter the element: ")
                    per3 = float(input("Percentage of this element: "))
                    mol1 = per1/float(df_atmass[df_atmass["Symbol"] == ele1].index.values)
                    mol2 = per2/float(df_atmass[df_atmass["Symbol"] == ele2].index.values)
                    mol3 = per3/float(df_atmass[df_atmass["Symbol"] == ele3].index.values)

                    if mol1<mol2 and mol1<mol3:
                        Mol1 = round(mol1/mol1)
                        Mol2 = round(mol2/mol1)
                        Mol3 = round(mol3/mol1)
                        print("The empirical formula is",ele1,ele2,Mol2,ele3,Mol3)

                    elif mol2<mol1 and mol2<mol3:
                        Mol1 = round(mol1/mol2)
                        Mol2 = 1
                        Mol3 = round(mol3/mol2)
                        print("The empirical formula is",ele1,Mol1,ele2,ele3,Mol3)

                    else:
                        Mol1 = round(mol1/mol3)
                        Mol2 = round(mol2/mol3)
                        Mol3 = 1
                        print("The empirical formula is",ele1,Mol1,ele2,Mol2,ele3)

                    mass_emp = (float(df_atmass[df_atmass["Symbol"] == ele1].index.values*Mol1) + float(df_atmass[df_atmass["Symbol"] == ele2].index.values*Mol2) + float(df_atmass[df_atmass["Symbol"] == ele3].index.values*Mol3))
                    emp_form = ele1,Mol1,ele2,Mol2,ele3,Mol3

                else:

                  ele1 = input("Enter the element: ")
                  per1 = float(input("Percentage of this element: "))
                  ele2 = input("Enter the element: ")
                  per2 = float(input("Percentage of this element: "))
                  mol1 = per1/float(df_atmass[df_atmass["Symbol"] == ele1].index.values)
                  mol2 = per2/float(df_atmass[df_atmass["Symbol"] == ele2].index.values)
                  if mol1<mol2:
                      Mol2 = round(mol2/mol1)
                      Mol1 = 1
                      print("The emperical formula is", ele1,ele2,Mol2)
                  else:
                      Mol1 = round(mol1/mol2)
                      Mol2 = 1
                      print("The emperical formula is",ele1,Mol1,ele2)

                  mass_emp = float((df_atmass[df_atmass["Symbol"] == df_atmass[df_atmass["Symbol"] == ele1].index.values].index.values*Mol1) + (df_atmass[df_atmass["Symbol"] == ele2].index.values*Mol2))
                  emp_form = ele1,Mol1,ele2,Mol2

                giv_mass = float(input("Enter given mass of compound: "))
                ratio = giv_mass/mass_emp
                print("The molecular formula of the compound is ",emp_form,ratio)

    elif chem_chp == "Atomic Structure" or chem_chp == "as":
      h = 6.626*(10**-34) 
      c = 3*(10**8)
      Na = 6.022*(10**23)
      Me = 9.11*(10**-31)
      Mp = 1.67*(10**-27)
      Mn = 1.67*(10**-27)
      pi = 3.14

      Help_atm = input("What do you need help with? Mass number (mn) , Wavelength , Frequency , Energy of photons (ep) , No. of photons emitted (npe) , KE of electron (ke) , Frequency of raditations emitted (fre) , Angular momentum of electron (ame) , Energy of sinlge electron species (esep) , Radius of single electron species (rsep) , Wavelength using de Broglie's equation (wdb), Mass using de Broglie's equation (mdb), Uncertainty in measurement (um) , Orbitals: ")
      if Help_atm == "Mass number" or Help_atm == "mass number" or Help_atm == "mn":
        print("Mass number is the sum of number of neutrons and number of protons")
        Massno = input('Enter the element of which you wish to find mass number:')
        mass_number = int(df_massno[df_massno["Symbol"] == Massno].index.values)
        print("Mass number is", mass_number)

      elif Help_atm == "Wavelength" or Help_atm == "wavelength":
        print("Wavelength w = c/v   where c = speed of electromagnetic radiation in vacuum, v = frequency")
        frequency = float(input("Enter frequency(in Hz): "))
        Wavelength = c/frequency
        print("Wavelength is", Wavelength,"m")
                
      elif Help_atm == "frequency" or Help_atm == "Frequency":
        print("Frequency v = c/w    where c = speed of electromagnetic radiation in vacuum, w = wavelength.")
        w = float(input("Enter wavelength(in nm)"))
        frequency = c/(w*(10**-9))
        print("Frequency is", frequency,"Hz")

      elif Help_atm == "Energy of photon" or Help_atm == "energy of photon" or Help_atm == "ep":
        print("Energy E = hv  where h = Planck'constant, v = frequency")
        v = float(input("Enter frequency(in Hz): "))
        E = h*v
        print("Energy of 1 photon is", E,"J")
        print("Energy of 1 mole of photons is", Na*E)
              
      elif Help_atm == "No. of photons emitted" or Help_atm == "no. of photons emitted" or Help_atm == "npe":
        print("No. of photons emitted = Power of bulb/Energy of photon")
        P = float(input("Enter power of bulb(in Watt): "))
        print("Energy of photon = h*v = (h*c)/w      where h = planck's constant, c =speed of electromagnetic radiation in vacuum, w = wavelength. ")
        given = input("Is frequency given or wavelength given?")
        if given == "frequency" or given == "Frequency":
          v = float(input("Enter frequency(in Hz): "))
          E = h*v
          print("Energy of photon is", E,"J")
          NPE = P/E
          print("No. of photons emitted is", NPE)
        else:
          w = float(input("Enter wavelength: "))
          E = (h*c)/(w*(10**-9))
          print("Energy of photon is", E)
          NPE = P/E
          print("No. of photons emitted is", NPE)

      elif Help_atm == "KE of electron" or Help_atm == "ke":
        print("KE of electron = mass of electron/ (frequency**2)  = h(v-vo)  where h = Planck's constant, v = frequency, vo = threshold frequency")
        v = float(input("Enter frequency: "))
        vo = float(input("Enter threshold frequency of metal: "))
        KE =h*(v-vo)
        print("Kinetic energy of electron is", KE)

      elif Help_atm == "Frequency of radiation emitted" or Help_atm == "frequency of radiation emitted" or Help_atm == "fre":
        print("Frequency of radiation emitted = difference in energy/Planck's constant ")
        fe = float(input("Enter final energy: "))
        ie = float(input("Enter initial energy: "))
        energy_diff = fe-ie
        fore = energy_diff/h
        print("Frequency of radiation emitted or absorbed is", fore)

      elif Help_atm == "Angular momentum of electron" or Help_atm == "angular momentum of electron" or Help_atm == "ame":
        print("Angular momentum of momentum = nh/2pi where n is the principal quantum number")
        n = int(input("Enter principal quantum number: "))
        AM = (n*h)/(2*pi)
        print("Angular momentum of electron is", AM)

      elif Help_atm == "Energy of single electron species" or Help_atm == "energy of single electron species" or Help_atm == "esep":
        print("Energies are given by this expression: En = -2.18 x 10**-18*(Z**2/n**2) where z = atomic number , n is principal quantum no.")
        Z = int(input("Enter atomic number of element: "))
        n = int(input("Enter the principal quantum number: "))    
        En = -2.18*(10**-18)*(Z**2/n**2)
        print("Energy of single electron species is", En)

      elif Help_atm == "Radius of single electron species" or Help_atm == "radius of single electron species " or Help_atm == "rsep":
        print("Expression: Rn = 52.9*((Z**2)/n")
        Z = int(input("Enter atomic number of element: "))
        n = int(input("Enter the principal quantum number: "))
        Rn = 52.9*((Z**2)/n)
        print("Radius of single electron species is", Rn)

      elif Help_atm == "Wavelength using de Broglie's equation" or Help_atm == "wavelength using de Broglie's equation" or Help_atm == "wdb":
        print("Expression: w = h/(m*v) = h/p where m is mass of particle, h is Planck's equation and v is frequency")
        m = float(input("Enter mass of particle: "))
        v = float(input("Enter frequency of particle: "))
        w = h/(m*v)
        print("Wavelength of particle is", w)

      elif Help_atm == "Mass using de Broglie's equation" or Help_atm == "mass using de Broglie's equation" or Help_atm == "mdb":
        print("Expression: m = h/(v*w) where w is wavelength of particle, v is frequency and h is Planck's constant" )
        v = float(input("Enter frequency of particle: "))
        w = float(input("Enter wavelength of particle: "))
        m = h/(v*w)
        print("Mass of particle is", m)

      elif Help_atm == "Uncertainty in measurement" or Help_atm == "uncertainty in measurement" or Help_atm == "um":
        print("According to Heisenberg's Uncertainty Principle: x*p = h/(4*pi*m) where x is uncertainty in postition and p is uncertainty in momentum")
        xorp = input("What do you want to find the uncertainty of? ")
        if xorp == "x":
          m = float(input("Enter mass of particle: "))
          p = float(input("Enter uncertainty in momentum: "))
          x = h/(4*pi*m*p)
          print("Uncertainty in position is", x)
        else:
          m = float(input("Enter mass of particle: "))
          x = float(input("Enter uncertainty in position: "))
          p = h/(4*pi*m*x)
          print("Uncertainty in momentum is", p)

      elif Help_atm == "Orbitals" or Help_atm == "orbitals":
        n = int(input("Enter principal quantum number: "))
        l = int(input("Enter Azimuthal quantum number: "))
        if l == 0:
          print("Orbital is {}s".format(n))
        elif l == 1:
          print("Orbital is {}p".format(n))
        elif l == 2:
          print("Orbital is {}d".format(n))
        elif l == 3:
          print("Orbital is {}f".format(n))
  
  else:
    print("Please enter valid subject.")

  quiz = input("Would you like to take a small test based on what you have learnt? (y/n)")
  if quiz == "y" or quiz == "yes" or quiz == "Y":
    sub = input("What subject do you want to take the quiz on? (P/C/M)")
    if sub == "M" or sub == "m" or sub == "math":
      import random
      chp = input("What Math chapter would you like to take a test for: Progressions (pr) or Straight lines(sl): ")
      #IDHAR PROGRESSIONS
      if chp == "Progressions" or chp == "progressions" or chp == "pr":
          num = random.randint(1,2)
          if num == 1:
              print("Q1) The 4 arithmetic means between 3 and 23 are: ")
              print("A) 5,9,11,13")
              print("B) 7,11,15,19")
              print("C) 5,11,15,22")
              print("D) 7,15,19,21")
              ans = input("Enter correct option: ")
              if ans == "B":
                print("Correct")
              else:
                print("Incorrect")
                  
              print()
              print("Q2) The GM of the numbers 3,(3^2),(3^3),...(3^n) is: ")
              print("A) 3^(2/n)")
              print("B) 3^((n+1)/2)")
              print("C) 3^(n/2)")
              print("D) 3^((n-1)/2)")
              ans = input("Enter correct option: ")
              if ans == "B":
                print("Correct")
              else:
                print("Incorrect")
                        
          else:
              print("Q1) The nth term of the series 3+10+17+... and 63+65+67+... are equal, then the value of n is?")
              print("A) 11")
              print("B) 12")
              print("C) 13")
              print("D) 15")
              ans = input("Enter correct option: ")
              if ans == "C":
                print("Correct")
              else:
                print("Incorrect")
              print()            
              print("Q2) The sum of few terms of any GP is 728, if common ratio is 3 and last term is 486, then first term of series will be?")
              print("A) 2")
              print("B) 1")
              print("C) 3")
              print("D) 4")
              ans = input("Enter correct option: ")
              if ans == "A":
                print("Correct")
              else:
                print("Incorrect")

      #IDHAR SE STRAIGHT LINES
      elif chp == "Straight lines" or chp == "sl" or chp == "straight lines":
          print("Q1) The equation of the line perpenicular to the line x/a - y/b = 1 and passing through the point at which it cuts x axis, is?")
          print("A) x/a + y/b + a/b = 0")
          print("B) x/b + y/a = b/a")
          print("C) x/b + y/a = 0")
          print("D) x/b + y/a = a/b")
          ans = input("Enter correct option: ")
          if ans == "A":
            print("Correct")
          else:
            print("Incorrect")
              
          print("Q2) Find the distance of the point (1,-1) from the line 12(x+6) = 5(y-2).")
          print("A) 4units")
          print("B) 8units")
          print("C) 6units")
          print("D) 5units")
          ans = input("Enter correct option: ")
          if ans == "D":
            print("Correct")
          else:
            print("Incorrect")

      else:
        print("Enter valid chapter")
        


    elif sub == "P" or sub == "p" or sub == "physics":
      chp = input("What physics chapter would you like to take the quiz for: Projectile Motion(pm) or Circular Motion(cm)?")
      if chp == "Projectile Motion" or chp == "pm":
        import random
        from PIL import Image
        num = random.randint(1,2)

        if num == 1:
            
            print('''Question 1. A particle is projected at an angle 37 deg with the incline plane in upward direction with speed 10 m/s. The angle of inclination of plane is 53 deg. Then the maximum distance from the incline plane
            attained by the particle will be:
            A)3m
            B)4m
            C)5m
            D)0m''')
            ans1 = input('Enter answer:')
            if ans1 == 'A':
              print("Good job! That's the correct answer!")
            else:
              print('''That answer is incorrect.''')
                
            print('''Question 2. It was calculated that a shell when fired from a gun with a certain velocity and at an angle of elevation 5pi/36 rad should strike a given target in same horizontal plane. In actual practice, it was found that a hill just prevented the trajectory. At what angle of elevation should the gun be
                  fired to hit the target:
                  A)5pi/36 rad
                  B)11pi/36 rad
                  C)7pi/36 rad
                  D)13pi/36 rad''')
            ans2 = input('Enter answer:')
            if ans2 == 'D':
              print("Good job that's the correct answer")
            else:              
              print("Incorrect")
                
            
        else:
            
            print('''Question 1. A point mass is projected, making an acute angle with the horizontal. If the angle between velocity vector and acceleration vector g is theta at any time
                  t during the motion, then theta is given by:
                  A)0 < theta < 90
                  B)theta = 90
                  C)theta < 90
                  D)0 < theta < 180''')
            ans3 = input("Enter answer:")
            if ans3 == 'D':
              print("Good job! That's the correct answer.")
            else:
              print("Incorrect")

            print('''Question 2. What is the maximum speed of oblique projectile from the ground in the vertical plane passing through a point (30m,40m) and projection
                  point is taken as the origin (g = 10 m/s^2):
                  A)30 m/s
                  B)20m/s
                  C)10root5 m/s
                  D)50 m/s''')
            ans4 = input("Enter answer:")
            if ans4 == "A":
              print("Good job! That answer's correct!")
            else:
              print("Incorrect")
      
      else:
        import random
        from PIL import Image
        num = random.randint(1,2)
        if num == 1:
          print('''Question 1. The maximum velocity with which a car driver must traverse a flat curve of radius 150m, coeff of friction 0.6 to avoid skidding
          A)60 m/s
          B)30 m/s
          C)15 m/s
          D)25 m/s''' )
          ans5 = input("Enter your answer:")
          if ans5 == "B":          
            print("Good job! That's the correct answer!")
          else:            
            print("Incorrect")
          print('''Question 2. A wheel is at rest. Its angular velocity increases uniformly and becomes 80 rad/s after 5 sec. Total angular displacement is:
                A)800 rad
                B)400 rad
                C)200 rad
                D)100 rad''')
          ans6 = input("Enter your answer:")
          if ans6 == 'C':
            print("Good job! That's the correct answer!")
          else:
            print("Incorrect")
          
        else:
          print('''Question 1. A particle moves along a circle of radius 20/pi m with constant tangential acceleration. If the speed of particle is 80 m/s at the end of the second revolution after the motion has begun, find tangential acceleration:
                A)160pi m/s^2
                B)40pi m/s^2
                C)40 m/s^2
                D)640pi m/s^2''')
          ans7 = input("Enter your answer:")
          if ans7 == "C":
            print("Good job! That's the correct answer!")
          else:
            print("Incorrect")
          print('''Question 2. A bucket's whirled in a vertical circle with a string. The water in bucket doesn't fall even when bucket's inverted at top of its path. In this position:
                A)mg = mv^2/r
                B)mg is greater than mv^2/r
                C)mg is not greater than mv^2/r
                D)mg is not less than mv^2/r''')
          ans8 = input("Enter your answer:")
          if ans8 == "C":
            print("Good job! That's the correct answer!")
          else:
            print("Incorrect")
      
    elif sub == "C" or sub == "c" or sub == "chemistry":
      import random
      chp = input("What Chemistry chapter would you like to take the quiz for: Mole Concept (mc) or Atomic Structure (as)?")
      if chp == "mc" or chp == "Mole Concept" or "mole concept":
          num = random.randint(1,2)
          if num == 1:
              print("Q1) Calculate the molarity of NaOH in the solution prepared by dissolving its 4 gms in enough water to form 250mL of the solution.")
              print("A) 0.3M")
              print("B) 4M")
              print("C) 0.4M")
              print("D) 3M")
              ans = input("Enter correct option: ")
              if ans == "C":
                  print("Correct")
              else:
                print("Incorrect")

              print()
                                  
              print("Q2) An organic compound contains 49.3% Carbon, 6.84% Hydrogen and its vapour density is 73. Molecular formula of compound is:")
              print("A) C3H5O2")
              print("B) C6H10O4")
              print("C) C3H10O2")
              print("D) C4H10O4")
              ans = input("Enter correct option: ")
              if ans == "B":               
                print("Correct")
              else:
                print("Incorrect")
                  

          else:
              print("Q1) The mole fraction of NaCl in a soltuion containing 1 mole of NaCl in 1000g of Water is: ")
              print("A) 0.0177")
              print("B) 0.001")
              print("C) 0.5")
              print("D) 1.5")
              ans = input("Enter correct option: ")
              if ans == "A":
                print("Correct")
              else:
                print("Incorrect")
              print()        
              print("Q2) A sample of clay was partially dried and then it contained 50% silica and 7% water. The original clay contained 12% water. Find the % of silica in the original sample.")
              print("A) 52.7%")
              print("B) 50%")
              print("C) 43%")
              print("D) 47.3%")
              ans = input("Enter correct option: ")
              if ans == "D":               
                print("Correct")
              else:
                print("Incorrect")
              
          
      elif chp == "as" or chp == "atomic structure" or chp == "Atomic Structure":
          num = random.randint(1,2)
          if num == 1:
              print("Q1) The energy of hydrogen atom in its ground state is -13.6eV. The energy of the level corresponding to n=5 is:")
              print("A) -0.54eV")
              print("B) -5.40eV")
              print("C) -0.85eV")
              print("D) -2.72eV")
              ans = input("Enter correct option: ")
              if ans == "A":
                print("Correct")
              else:
                print("Incorrect")
              print()        
              print("Q2) de-Broglie wavelength of electron in second orbit of Li2+ ion will be equal to de-Broglie of wavelength of electron in:")
              print("A) n=3 of H-atom")
              print("B) n=4 of C5+ atom")
              print("C) n=6 of Be3+ atom")
              print("D) n=3 of He+ atom")
              ans = input("Enter correct option: ")
              if ans == "B":
                print("Correct")
              else:
                print("Incorrect")

          else:
              print("Q1 The frequency of yellow light having wavelength 600nm is:")
              print("A) 5 x 10^14 Hz")
              print("B) 2.5 x 10^7 Hz")
              print("C) 5 x 10^7 Hz")
              print("D) 2.5 x 10^14 Hz")
              ans = input("Enter correct option: ")
              if ans == "A":            
                print("Correct")
              else:
                print("Incorrect")
              print()        
              print("Q2) The uncertainty in the moment of an electron is 1 x 10^-5 kgm/s. The uncertainty in its position will be: (h = 6.626 x 10^-34Js)")
              print("A) 1.05 x 10^-28 m")
              print("B) 1.05 x 10^-26 m")
              print("C) 5.27 x 10^-30 m")
              print("D) 5.25 x 10^-28 m")
              ans = input("Enter correct option: ")
              if ans == "C":
                print("Correct")
              else:
                print("Incorrect")
      
      else:
          print("Enter valid chapter")

      

    else:
      print("Enter valid subject")

  else:
    print("Happy learning!")

    


    
    
