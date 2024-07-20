# FITexel

## Overview
This project implements a simple spreadsheet evaluation system in C++. 
The core functionality includes support for various cell types (numeric, string, references), as well as the evaluation of expressions involving arithmetic and logical operations.

## Project Structure
- CPos
- Node
  - BinaryOperatorNode
    - AddNode  (operator +)
    - SubNode  (operator -)
    - MultNode (operator *)
    - DivNode  (operator /)
    - PowNode  (operator ^)
    - EqNode   (operator ==)
    - NeqNode  (operator !=)
    - LtNode   (operator <)
    - LeNode   (operator <=)
    - GtNode   (operator >)
    - GeNode   (operator >=)
  - ValueNode
  - ReferenceNode
  - NegNode    (operator unary -)
- CExprBuilder

## Classes and Functionality

### Cpos
Represents a cell position in the spreadsheet.

### Node
Abstract base class for all types of nodes in the expression tree. Includes methods for evaluating expressions, cloning nodes, and serializing data into a stream.

### BinaryOperatorNode
Base class for all binary operation nodes. Subclass of the Node class.

### ValueNode
Represents a cell containing a value (undefined, string, or numeric).

### ReferenceNode
Represents a cell reference, allowing for relative and absolute references.

### AddNode, SubNode, MultNode, DivNode, PowNode, EqNode, NeqNode, LtNode, LeNode, GtNode, GeNode
Concrete subclasses of BinaryOperatorNode, implementing specific arithmetic and logical operations.

### NegNode
Subclass of Node, implements unary negation operator.

### CExprBuilder
A helper class to parse expressions. When parseExpression function goes through an expression, it translates it from infix representation into postfix by invoking corresponding methods of CExprBuilder. 
For example, in case of expression ```1 + 2```, parsing it will result in the following methods call chain of CExprBuilder: ```valNumber (1), valNumber(2), opAdd ()```.   

This project is my university project, so the real implementation of parseExpression was hidden in a testing environment and to us it was provided as a statically linked library (libexpression_parser.a). 
As for CExprBuilder class, it is an abstract class and the real implementation is provided via it's subclass: ExpressionBuilder.

## CSpreadsheet Class and It's Functionality
This class is the spreadsheet processor itself. It provides the following interface:

- **capabilities**: method used by the testing environment, messages what additional features were implemented (in my case, cyclic dependencies detection nd parsing process optimization).

- **Default constructor**: constructs an emptu spreadsheet.

- **Copy constructor**: creates a deep copy of a spreadsheet.

- **Move constructor**: creates a spreadsheet copy, while also resets to the default state the source spreadsheet instance.

- **Move and Copy equivalence operators**: implements copying and moving by copy-and-swap idiom.

- **load**: method loads a spreadsheet from a stream.

- **save**: method save a spreadsheet into a stream.

- **setCell**: method sets the cell value inside a spreadsheet. To do so, it requires a position and the content parameters.

- **getValue**: returns an evaluated value of a cell at the provided position. In case of cyclic dependency, throws a corresponding exception.

- **copyRect** copies a rectangle from a source position into a destination position. Rectangle dimensions are provided as w(width) and h(height) parameters.

## Summary
This project has been done in my university and due to tight deadlines, it's been done in a bit of a hurry, so it lacks clarifying comments. 
Also, due to a testing system, it required to put all the classes in one file, which is also not beneficial for the code quality, so in the future some minor refactoring has to be done and the code should be separated into different files.























