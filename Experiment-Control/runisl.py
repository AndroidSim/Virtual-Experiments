#!/usr/bin/python

import os
import javabridge

##
##  without running "ant jar" to create the dist/isl.jar
##
javaclasspath = [os.path.realpath(os.path.join(os.path.dirname(__file__),'..','build','classes'))] + \
                [os.path.realpath(os.path.join(os.path.dirname(__file__),'..','cfg'))] + \
                [os.path.realpath(os.path.join(os.path.dirname(__file__),'..','lib', name + '.jar'))
                   for name in ['hamcrest-core-1.1', 'jmf', 'colt-1.2.0', 'iText-2.1.5', 'junit4-4.8.2', 
                                'commons-math-2.2', 'jcommon-1.0.16', 'logback-classic-0.9.28', 
                                'mason17-with-src', 'ecj19', 'jfreechart-1.0.13', 'logback-core-0.9.28', 
                                'slf4j-api-1.6.1']] + \
                javabridge.JARS

##
## using dist/isl.jar
##
#javaclasspath = [os.path.realpath(os.path.join(os.path.dirname(__file__),'..','dist', 'isl.jar'))] + \
#                [os.path.realpath(os.path.join(os.path.dirname(__file__),'..','cfg'))] + \
#                javabridge.JARS

#print os.pathsep.join(javaclasspath)

javabridge.start_vm(class_path=[os.pathsep.join(javaclasspath)], run_headless=True)

try:
    print javabridge.get_env().get_version()
    print javabridge.run_script('java.lang.System.getProperty("java.class.path")')

    ###
    ## rhino invocation of isl.Main.main(String[] args)
    ### 
    #print javabridge.run_script('Packages.isl.Main.main([])');


    ###
    ## high invocation of isl.Main.main(String[] args)
    ###
    #javabridge.static_call("isl/Main", "main", "([Ljava/lang/String;)V", [])

    ###
    ## low API invocation of isl.Main.main(String[] args)
    ###
    main_class = javabridge.get_env().find_class("isl/Main");
    main_method_id = javabridge.get_env().get_static_method_id(main_class,'main','([Ljava/lang/String;)V')

    ###### ways to construct the String[] argument

    ## rhino construction of String[]
    #main_args = javabridge.run_script("[]")  # cheat by creating the array in rhino

    ## high construction of String[]
    #main_args = javabridge.make_list()

    ### low construction of String[]
    string_class = javabridge.get_env().find_class("java/lang/String")
    main_args = javabridge.get_env().make_object_array(0,string_class)


    ### now that we have the String[] arg in a Java Object form we can invoke!
    javabridge.get_env().call_static_method(main_class,main_method_id,main_args)



    ###
    ## other Examples
    ###
    
    ### Use the high and low level apis to create an Integer and a MyInt
    #print('Using api: ')
    #mic = javabridge.get_env().find_class("isl/util/MyInt")
    #mi = javabridge.make_instance("isl/util/MyInt", "(J)V", 12)
    #print javabridge.call(mi,"doubleValue","()D")
    #
    #i = javabridge.make_instance("java/lang/Integer","(I)V", 12)
    #print javabridge.call(i,"toString","()Ljava/lang/String;")

    ### Use Rhino (JavaScript) to create a MyInt
    #print('Using Rhino: ')
    #print javabridge.run_script('mi = new Packages.isl.util.MyInt(12);mi.doubleValue();')

finally:
    javabridge.kill_vm()
