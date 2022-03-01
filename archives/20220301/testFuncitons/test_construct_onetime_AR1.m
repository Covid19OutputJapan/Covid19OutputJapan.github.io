classdef test_construct_onetime_AR1 < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        function addPath(testCase)
            addpath('..')
        end
    end
    
    methods(Test)
        function example1(testCase)
            actSolution =  construct_onetime_AR1(0.5, 0.1, 4) ;
            expSolution = [1.1; 1.05; 1.025; 1.0125];
            testCase.verifyEqual(actSolution,expSolution)
        end
        function example2(testCase)
            actSolution =  construct_onetime_AR1(0.5, -0.1, 4) ;
            expSolution = [0.9; 0.95; 0.975; 0.9875];
            testCase.verifyEqual(actSolution,expSolution)
        end
        function inputPersistenceErrorNegativeReal(testCase)
            testCase.verifyError(@()construct_onetime_AR1(-1, -0.1, 4), "MATLAB:NotPositiveReal")
        end
        function inputPersistenceErrorNonNumeric(testCase)
            testCase.verifyError(@()construct_onetime_AR1('1', -0.1, 4), "MATLAB:NotPositiveReal")
        end
        function inputShockSizeErrorNonNumeric(testCase)
            testCase.verifyError(@()construct_onetime_AR1(0.5, '10', 4), "MATLAB:NotReal")
        end
        function inputShockSizeErrorImmaginary(testCase)
            testCase.verifyError(@()construct_onetime_AR1(0.5, 0+0.1*1i, 4), "MATLAB:NotReal")
        end
        function inputLengthErrorNegativeInteger(testCase)
            testCase.verifyError(@()construct_onetime_AR1(0.5, -0.9, -5),"MATLAB:NotPositiveInteger")
        end
        function inputLengthErrorNonInteger(testCase)
            testCase.verifyError(@()construct_onetime_AR1(0.5,  -0.9, 9.99), "MATLAB:NotPositiveInteger")
        end
        function inputLengthErrorNonNumeric(testCase)
            testCase.verifyError(@()construct_onetime_AR1(0.5, -0.9, '-5'),"MATLAB:NotPositiveInteger")
        end
    end
end