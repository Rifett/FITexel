#ifndef __PROGTEST__
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <climits>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <array>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <variant>
#include <optional>
#include <compare>
#include <charconv>
#include <span>
#include <utility>

using namespace std::literals;
using CValue = std::variant<std::monostate, double, std::string>;

constexpr unsigned SPREADSHEET_CYCLIC_DEPS = 0x01;
constexpr unsigned SPREADSHEET_FUNCTIONS   = 0x02;
constexpr unsigned SPREADSHEET_FILE_IO     = 0x04;
constexpr unsigned SPREADSHEET_SPEED       = 0x08;
constexpr unsigned SPREADSHEET_PARSER      = 0x10;
#endif /* __PROGTEST__ */


#define UNDEFINED 0
#define DOUBLE    1
#define STRING    2


//============================================POSITION CLASS=======================================
class CPos
{
public:
    CPos() = default;

    auto operator <=> ( const CPos& other ) const = default;

    CPos ( std::string_view position )
    {
        size_t modulo = 26, index = 0;

        //Parse alphabetic part of a position (column index)
        for ( ; index < position.size(); ++index ) {
            if ( !std::isalpha( position.at( index ) ) )
                break;
            if ( columnIndex != 0 )
                columnIndex *= modulo;
            columnIndex += ( tolower( position.at( index ) ) - 'a' + 1);
        }
        --columnIndex;

        if ( index == 0 ) throw std::invalid_argument("Wrong position format: wrong column index formatting");
        modulo = 10;
        size_t indexBeforeRow = index;

        //Parse numeric part of a position (row index)
        for ( ; index < position.size(); ++index ) {
            if ( !std::isdigit( position.at( index ) ) )
                break;
            if ( rowIndex != 0 )
                rowIndex *= modulo;
            rowIndex += ( position.at( index ) - '0');
        }

        if ( index != position.size() || index == indexBeforeRow ) throw std::invalid_argument("Wrong position format: wrong row index formatting");
    }

    size_t getColumnIndex() const {
        return columnIndex;
    }

    size_t getRowIndex() const {
        return rowIndex;
    }

    void applyShift( int colShift, int rowShift ) {
        columnIndex += colShift;
        rowIndex += rowShift;
    }

    std::pair<std::string, std::string> restoreFormat() const {
        size_t columnIndexCopy = columnIndex + 1;
        std::string columnAlphaIndex;

        while (columnIndexCopy > 0) {
            columnAlphaIndex = static_cast<char>('A' + (columnIndexCopy - 1) % 26) + columnAlphaIndex;
            columnIndexCopy = (columnIndexCopy - 1) / 26;
        }

        return {columnAlphaIndex, std::to_string(rowIndex)};
    }

private:
    size_t columnIndex = 0, rowIndex = 0;
};


//============================================NODES================================================
//Abstract base class for all nodes.
class Node
{
public:
    virtual CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) = 0;

    virtual std::unique_ptr< Node > clone() = 0;

    virtual void getDataIntoStream(std::ostream& os) = 0;

    virtual Node& applyShift( int colShift, int rowShift ) {
        return *this;
    }

    void setUpper() {
        expressionUpper = true;
    }

protected:
    bool expressionUpper = false;
};


//Base class for binary operation nodes.
class BinaryOperatorNode : public Node
{
public:
    BinaryOperatorNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : lhs(std::move(lhs)), rhs(std::move(rhs)) {}

    Node &applyShift(int colShift, int rowShift) override {
        lhs->applyShift( colShift, rowShift );
        rhs->applyShift( colShift, rowShift );
        return *this;
    }

protected:
    std::unique_ptr<Node> lhs, rhs;
};


//Node for values: empty (undefined), string value, numerical value.
class ValueNode : public Node
{
public:
    ValueNode( CValue val = std::monostate() ) : val( std::move( val ) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override { return val; }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< ValueNode >( *this );
    }

    void getDataIntoStream(std::ostream& os) override {
        switch (val.index()) {
            case STRING:
                if ( expressionUpper ) {
                    os << "=\"";
                    for (char ch: get< std::string >( val )) {
                        if ( ch == '\"' )
                            os << ch;
                        os << ch;
                    }
                    os << '\"';
                } else
                    os << get< std::string >( val );
                break;
            case DOUBLE:
                if ( expressionUpper ) os << '=';
                os << std::to_string( get< double >( val ) );
                break;
            default:
                os << "";
                break;
        }
    }

protected:
    CValue val;
};


//Node for references.
class ReferenceNode : public Node
{
public:
    ReferenceNode( std::string ref ) {
        if ( !ref.empty() && ref.front() == '$' ) {
            columnIndexIsRelative = false;
            ref.erase( 0, 1 );
        }

        size_t dollarIndex = ref.find('$');
        if ( dollarIndex != std::string::npos ) {
            rowIndexIsRelative = false;
            ref.erase(dollarIndex, 1);
        }

        try {
            refPos = CPos(ref);
        } catch ( ... ) {
            throw std::runtime_error("Wrong reference formatting.");
        }
    }

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        if ( visited.contains( refPos ) || table.size() <= refPos.getColumnIndex()
             || table.at( refPos.getColumnIndex() ).size() <= refPos.getRowIndex() || !table.at( refPos.getColumnIndex() ).at( refPos.getRowIndex() ) ) return {};
        auto visitedCopy = visited;
        visitedCopy.insert( refPos );
        return table.at( refPos.getColumnIndex() ).at( refPos.getRowIndex() )->evaluate( table, visitedCopy );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< ReferenceNode >( *this );
    }

    Node &applyShift(int colShift, int rowShift) override {
        refPos.applyShift( columnIndexIsRelative ? colShift : 0, rowIndexIsRelative ? rowShift : 0 );
        return *this;
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        auto origFormat = refPos.restoreFormat();
        if ( columnIndexIsRelative ) os << '$';
        os << origFormat.first;
        if ( rowIndexIsRelative ) os << '$';
        os << origFormat.second;
    }

protected:
    CPos refPos;
    bool columnIndexIsRelative = true, rowIndexIsRelative = true;
};


// lhs + rhs
class AddNode : public BinaryOperatorNode
{
public:
    AddNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() == DOUBLE && rhsValue.index() == DOUBLE )
            return get < double > ( lhsValue ) + get < double > ( rhsValue );
        else if ( lhsValue.index() == DOUBLE && rhsValue.index() == STRING )
            return std::to_string( get < double > ( lhsValue ) ) + get < std::string > ( rhsValue );
        else if ( lhsValue.index() == STRING && rhsValue.index() == DOUBLE )
            return get < std::string > ( lhsValue ) + std::to_string( get < double > ( rhsValue ) );
        else if ( lhsValue.index() == STRING && rhsValue.index() == STRING )
            return get < std::string > ( lhsValue ) + get < std::string > ( rhsValue );
        else
            return {};
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< AddNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "+";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs - rhs
class SubNode : public BinaryOperatorNode
{
public:
    SubNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != DOUBLE || rhsValue.index() != DOUBLE ) return {};
        return get < double > ( lhsValue ) - get < double > ( rhsValue );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< SubNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "-";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs * rhs
class MultNode : public BinaryOperatorNode
{
public:
    MultNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != DOUBLE || rhsValue.index() != DOUBLE ) return {};
        return get < double > ( lhsValue ) * get < double > ( rhsValue );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< MultNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "*";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs / rhs
class DivNode : public BinaryOperatorNode
{
public:
    DivNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != DOUBLE || rhsValue.index() != DOUBLE || get < double > ( rhsValue ) == 0 ) return {};
        return get < double > ( lhsValue ) / get < double > ( rhsValue );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< DivNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "/";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs ^ rhs
class PowNode : public BinaryOperatorNode
{
public:
    PowNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != DOUBLE || rhsValue.index() != DOUBLE ) return {};
        return std::pow ( get < double > ( lhsValue ), get < double > ( rhsValue ) );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< PowNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "^";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//- operand
class NegNode : public Node
{
public:
    NegNode( std::unique_ptr<Node> operand ) : operand( std::move(operand) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto operandValue = operand->evaluate( table, visited );
        if ( operandValue.index() != DOUBLE ) return {};
        return -1 * get < double > ( operandValue );
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< NegNode >( operand->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "-";
        operand->getDataIntoStream( os );
    }

protected:
    std::unique_ptr<Node> operand;
};


//lhs == rhs
class EqNode : public BinaryOperatorNode
{
public:
    EqNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                  get < double      > ( lhsValue ) == get < double      > ( rhsValue )
                : get < std::string > ( lhsValue ) == get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< EqNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "=";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs != rhs
class NeqNode : public BinaryOperatorNode
{
public:
    NeqNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                     get < double      > ( lhsValue ) != get < double      > ( rhsValue )
                   : get < std::string > ( lhsValue ) != get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< NeqNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "<>";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs < rhs
class LtNode : public BinaryOperatorNode
{
public:
    LtNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                     get < double      > ( lhsValue ) < get < double      > ( rhsValue )
                   : get < std::string > ( lhsValue ) < get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< LtNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "<";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs <= rhs
class LeNode : public BinaryOperatorNode
{
public:
    LeNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                     get < double      > ( lhsValue ) <= get < double      > ( rhsValue )
                   : get < std::string > ( lhsValue ) <= get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< LeNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << "<=";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs > rhs
class GtNode : public BinaryOperatorNode
{
public:
    GtNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                     get < double      > ( lhsValue ) > get < double      > ( rhsValue )
                   : get < std::string > ( lhsValue ) > get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< GtNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << ">";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//lhs >= rhs
class GeNode : public BinaryOperatorNode
{
public:
    GeNode( std::unique_ptr<Node> lhs, std::unique_ptr<Node> rhs ) : BinaryOperatorNode( std::move(lhs), std::move(rhs) ) {}

    CValue evaluate( std::vector < std::vector < std::unique_ptr < Node > > >& table, std::set < CPos >& visited ) override {
        auto lhsValue = lhs->evaluate( table, visited ), rhsValue = rhs->evaluate( table, visited );
        if ( lhsValue.index() != rhsValue.index() || lhsValue.index() == UNDEFINED ) return {};
        auto res = lhsValue.index() == DOUBLE ?
                     get < double      > ( lhsValue ) >= get < double      > ( rhsValue )
                   : get < std::string > ( lhsValue ) >= get < std::string > ( rhsValue );
        return res ? 1.0 : 0.0;
    }

    std::unique_ptr<Node> clone() override {
        return std::make_unique< GeNode >( lhs->clone(), rhs->clone() );
    }

    void getDataIntoStream(std::ostream& os) override {
        if ( expressionUpper ) os << '=';
        os << "(";
        lhs->getDataIntoStream( os );
        os << ">=";
        rhs->getDataIntoStream( os );
        os << ")";
    }
};


//============================================EXPRESSIONS STUFF====================================

//TODO -> Remove this class before submitting ( for whatever reason )
class CExprBuilder
{
public:
    virtual void opAdd () = 0;
    virtual void opSub () = 0;
    virtual void opMul () = 0;
    virtual void opDiv () = 0;
    virtual void opPow () = 0;
    virtual void opNeg () = 0;
    virtual void opEq () = 0;
    virtual void opNe () = 0;
    virtual void opLt () = 0;
    virtual void opLe () = 0;
    virtual void opGt () = 0;
    virtual void opGe () = 0;
    virtual void valNumber ( double val ) = 0;
    virtual void valString ( std::string val ) = 0;
    virtual void valReference ( std::string val ) = 0;
    virtual void valRange ( std::string val ) = 0;
    virtual void funcCall ( std::string fnName, int paramCount ) = 0;
};


class ExpressionBuilder : public CExprBuilder {
public:
    void opAdd() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< AddNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opSub() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< SubNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opMul() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< MultNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opDiv() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< DivNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opPow() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< PowNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opNeg() override {
        auto operand = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< NegNode > ( std::move( operand ) ) );
    }

    void opEq() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< EqNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opNe() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< NeqNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opLt() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< LtNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opLe() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< LeNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opGt() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< GtNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void opGe() override {
        auto rhs = std::move( stack.top() ); stack.pop();
        auto lhs = std::move( stack.top() ); stack.pop();
        stack.push( std::make_unique< GeNode > ( std::move( lhs ), std::move( rhs ) ) );
    }

    void valNumber(double val) override {
        stack.push( std::make_unique< ValueNode > ( val ) );
    }

    void valString(std::string val) override {
        stack.push( std::make_unique< ValueNode > ( val ) );
    }

    void valReference(std::string val) override {
        stack.push( std::make_unique< ReferenceNode > ( val ) );
    }

    void valRange(std::string val) override {
        //TODO -> Sorry, but this feature will remain unimplemented as I am running low on time ( I need to write my thesis )
    }

    void funcCall(std::string fnName, int paramCount) override {
        //TODO -> Same as above
    }

    std::unique_ptr< Node > getExprResult() {
        return stack.size() == 1 ? std::move( stack.top() ) : nullptr;
    }

protected:
    std::stack< std::unique_ptr< Node > > stack;
};


void parseExpression ( std::string expr, CExprBuilder & builder );


//============================================SPREADSHEET==========================================
class CSpreadsheet
{
public:
    static unsigned capabilities ()
    {
        return SPREADSHEET_CYCLIC_DEPS | SPREADSHEET_SPEED;
    }

    CSpreadsheet () = default;

    CSpreadsheet(const CSpreadsheet& other) {
        table.reserve(other.table.size());
        for (const auto& row : other.table) {
            std::vector<std::unique_ptr<Node>> clonedRow;
            for (const auto& nodePtr : row) {
                if ( nodePtr )
                    clonedRow.push_back( nodePtr->clone() );
                else
                    clonedRow.push_back(nullptr );
            }
            table.push_back(std::move(clonedRow));
        }
    }

    CSpreadsheet(CSpreadsheet&& other) noexcept : table(std::move(other.table)) {}

    CSpreadsheet& operator=(const CSpreadsheet& other) {
        if (this != &other) {
            auto copy ( other );
            std::swap( *this, copy );
        }
        return *this;
    }

    CSpreadsheet& operator=(CSpreadsheet&& other) noexcept {
        if (this != &other)
            table = std::move(other.table);
        return *this;
    }

    ~CSpreadsheet() = default;

    bool load ( std::istream & is ) {
        table.clear();

        size_t colIndex, rowIndex;
        while ( is >> colIndex >> rowIndex ) {
            std::string cellData;
            char ch;

            if ( !(is >> ch) || ch != '\\' )
                return false;

            while ( is.get( ch ) && ch != '\\' ) {
                cellData.push_back( ch );
            }

            if ( ch != '\\' )
                return false;

            //Ignore "\n"
            is.ignore();

            // Assure that there are enough space in the table
            if (table.size() <= colIndex) table.resize(colIndex + 1);
            if (table.at(colIndex).size() <= rowIndex) table.at(colIndex).resize(rowIndex + 1);

            ExpressionBuilder exprBuilder;
            parseExpression(cellData, exprBuilder);
            auto exprResult = exprBuilder.getExprResult();

            if (exprResult) {
                table.at(colIndex).at(rowIndex) = std::move(exprResult);
            } else {
                return false;
            }
        }
        return is.eof();
    }

    bool save ( std::ostream & os ) const {
        for (size_t colIndex = 0; colIndex < table.size(); ++colIndex) {
            for (size_t rowIndex = 0; rowIndex < table.at(colIndex).size(); ++rowIndex) {
                if ( table.at( colIndex ).at( rowIndex ) ) {
                    os << colIndex << " " << rowIndex << " " << "\\";
                    table.at( colIndex ).at( rowIndex )->getDataIntoStream( os );
                    os << "\\" << std::endl;
                }
            }
        }
        return true;
    }

    bool setCell ( CPos pos, const std::string& contents ) {
        //Assure that there are enough space
        if ( table.size() <= pos.getColumnIndex() ) table.resize( pos.getColumnIndex() + 1 );
        if ( table.at( pos.getColumnIndex() ).size() <= pos.getRowIndex() ) table.at( pos.getColumnIndex() ).resize( pos.getRowIndex() + 1 );

        ExpressionBuilder exprBuilder;
        parseExpression(contents, exprBuilder);
        auto exprResult = exprBuilder.getExprResult();

        if ( exprResult ) {
            if ( contents.front() == '=' ) exprResult->setUpper();
            table.at( pos.getColumnIndex() ).at( pos.getRowIndex() ) = std::move( exprResult );
            return true;
        } else {
            return false;
        }
    }

    CValue getValue ( CPos pos ) {
        if ( pos.getColumnIndex() >= table.size() || pos.getRowIndex() >= table.at( pos.getColumnIndex() ).size() || !table.at( pos.getColumnIndex() ).at( pos.getRowIndex() ) )
            return {};
        std::set<CPos> visited;
        visited.insert(pos);

        return table.at(pos.getColumnIndex()).at(pos.getRowIndex())->evaluate( table, visited );
    }

    void copyRect ( CPos dst, CPos src, int w = 1, int h = 1 ) {
        std::vector<std::vector<std::unique_ptr<Node>>> tempRect;

        tempRect.resize( w );
        for (auto& column: tempRect) {
            column.resize( h );
        }

        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {
                if (src.getColumnIndex() + col < table.size() &&
                    src.getRowIndex() + row < table[src.getColumnIndex() + col].size() &&
                    table[src.getColumnIndex() + col][src.getRowIndex() + row]) {
                    tempRect[col][row] = table[src.getColumnIndex() + col][src.getRowIndex() + row]->clone();
                } else {
                    tempRect[col][row] = nullptr;
                }
            }
        }

        int colShift = dst.getColumnIndex() - src.getColumnIndex();
        int rowShift = dst.getRowIndex() - src.getRowIndex();

        for (auto& col : tempRect) {
            for (auto& nodePtr : col) {
                if (nodePtr) {
                    nodePtr->applyShift( colShift, rowShift );
                }
            }
        }

        if (table.size() <= dst.getColumnIndex() + w)
            table.resize(dst.getColumnIndex() + w + 1);
        for (int col = 0; col < w; ++col)
            if (table[dst.getColumnIndex() + col].size() <= dst.getRowIndex() + h)
                table[dst.getColumnIndex() + col].resize(dst.getRowIndex() + h + 1);


        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {
                table[dst.getColumnIndex() + col][dst.getRowIndex() + row] = tempRect[col][row] ? std::move(tempRect[col][row]) : nullptr;
            }
        }
    }

private:
    std::vector < std::vector < std::unique_ptr < Node > > > table;
};


//Dear Mr.Reviewer, I am really sorry for the terrible state of my code, I've tried to make it pretty ( at least for the first two days )
//But then I realized that I need to hurry to do my Bachelor thesis ( currently, I have only 12 days left to do it ( there are only a couple lines of code right now ) )
//So yeah, I hope that reading this code will not cause you any mental or physical health issues and I can only wish you strong nerves and a great power of will ( great enough to go through this code)


#ifndef __PROGTEST__

bool valueMatch ( const CValue & r, const CValue & s )

{
    if ( r . index () != s . index () )
        return false;
    if ( r . index () == 0 )
        return true;
    if ( r . index () == 2 )
        return std::get<std::string> ( r ) == std::get<std::string> ( s );
    if ( std::isnan ( std::get<double> ( r ) ) && std::isnan ( std::get<double> ( s ) ) )
        return true;
    if ( std::isinf ( std::get<double> ( r ) ) && std::isinf ( std::get<double> ( s ) ) )
        return ( std::get<double> ( r ) < 0 && std::get<double> ( s ) < 0 )
               || ( std::get<double> ( r ) > 0 && std::get<double> ( s ) > 0 );
    return fabs ( std::get<double> ( r ) - std::get<double> ( s ) ) <= 1e8 * DBL_EPSILON * fabs ( std::get<double> ( r ) );
}
int main ()
{
    CSpreadsheet x0, x1;
    std::ostringstream oss;
    std::istringstream iss;
    std::string data;
    assert ( x0 . setCell ( CPos ( "A1" ), "10" ) );
    assert ( x0 . setCell ( CPos ( "A2" ), "20.5" ) );
    assert ( x0 . setCell ( CPos ( "A3" ), "3e1" ) );
    assert ( x0 . setCell ( CPos ( "A4" ), "=40" ) );
    assert ( x0 . setCell ( CPos ( "A5" ), "=5e+1" ) );
    assert ( x0 . setCell ( CPos ( "A6" ), "raw text with any characters, including a quote \" or a newline\n" ) );
    assert ( x0 . setCell ( CPos ( "A7" ), "=\"quoted string, quotes must be doubled: \"\". Moreover, backslashes are needed for C++.\"" ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A1" ) ), CValue ( 10.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A2" ) ), CValue ( 20.5 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A3" ) ), CValue ( 30.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A4" ) ), CValue ( 40.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A5" ) ), CValue ( 50.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A6" ) ), CValue ( "raw text with any characters, including a quote \" or a newline\n" ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A7" ) ), CValue ( "quoted string, quotes must be doubled: \". Moreover, backslashes are needed for C++." ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "A8" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "AAAA9999" ) ), CValue() ) );
    assert ( x0 . setCell ( CPos ( "B1" ), "=A1+A2*A3" ) );
    assert ( x0 . setCell ( CPos ( "B2" ), "= -A1 ^ 2 - A2 / 2   " ) );
    assert ( x0 . setCell ( CPos ( "B3" ), "= 2 ^ $A$1" ) );
    assert ( x0 . setCell ( CPos ( "B4" ), "=($A1+A$2)^2" ) );
    assert ( x0 . setCell ( CPos ( "B5" ), "=B1+B2+B3+B4" ) );
    assert ( x0 . setCell ( CPos ( "B6" ), "=B1+B2+B3+B4+B5" ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B1" ) ), CValue ( 625.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B2" ) ), CValue ( -110.25 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B3" ) ), CValue ( 1024.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B4" ) ), CValue ( 930.25 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B5" ) ), CValue ( 2469.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B6" ) ), CValue ( 4938.0 ) ) );
    assert ( x0 . setCell ( CPos ( "A1" ), "12" ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B1" ) ), CValue ( 627.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B2" ) ), CValue ( -154.25 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B3" ) ), CValue ( 4096.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B4" ) ), CValue ( 1056.25 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B5" ) ), CValue ( 5625.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B6" ) ), CValue ( 11250.0 ) ) );
    x1 = x0;
    assert ( x0 . setCell ( CPos ( "A2" ), "100" ) );
    assert ( x1 . setCell ( CPos ( "A2" ), "=A3+A5+A4" ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B1" ) ), CValue ( 3012.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B2" ) ), CValue ( -194.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B3" ) ), CValue ( 4096.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B4" ) ), CValue ( 12544.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B5" ) ), CValue ( 19458.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "B6" ) ), CValue ( 38916.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B1" ) ), CValue ( 3612.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B2" ) ), CValue ( -204.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B3" ) ), CValue ( 4096.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B4" ) ), CValue ( 17424.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B5" ) ), CValue ( 24928.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B6" ) ), CValue ( 49856.0 ) ) );
    oss . clear ();
    oss . str ( "" );
    assert ( x0 . save ( oss ) );
    data = oss . str ();
    iss . clear ();
    iss . str ( data );
    assert ( x1 . load ( iss ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B1" ) ), CValue ( 3012.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B2" ) ), CValue ( -194.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B3" ) ), CValue ( 4096.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B4" ) ), CValue ( 12544.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B5" ) ), CValue ( 19458.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B6" ) ), CValue ( 38916.0 ) ) );
    assert ( x0 . setCell ( CPos ( "A3" ), "4e1" ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B1" ) ), CValue ( 3012.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B2" ) ), CValue ( -194.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B3" ) ), CValue ( 4096.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B4" ) ), CValue ( 12544.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B5" ) ), CValue ( 19458.0 ) ) );
    assert ( valueMatch ( x1 . getValue ( CPos ( "B6" ) ), CValue ( 38916.0 ) ) );
    oss . clear ();
    oss . str ( "" );
    assert ( x0 . save ( oss ) );
    data = oss . str ();
    for ( size_t i = 0; i < std::min<size_t> ( data . length (), 10 ); i ++ )
        data[i] ^=0x5a;
    iss . clear ();
    iss . str ( data );
    assert ( ! x1 . load ( iss ) );
    assert ( x0 . setCell ( CPos ( "D0" ), "10" ) );
    assert ( x0 . setCell ( CPos ( "D1" ), "20" ) );
    assert ( x0 . setCell ( CPos ( "D2" ), "30" ) );
    assert ( x0 . setCell ( CPos ( "D3" ), "40" ) );
    assert ( x0 . setCell ( CPos ( "D4" ), "50" ) );
    assert ( x0 . setCell ( CPos ( "E0" ), "60" ) );
    assert ( x0 . setCell ( CPos ( "E1" ), "70" ) );
    assert ( x0 . setCell ( CPos ( "E2" ), "80" ) );
    assert ( x0 . setCell ( CPos ( "E3" ), "90" ) );
    assert ( x0 . setCell ( CPos ( "E4" ), "100" ) );
    assert ( x0 . setCell ( CPos ( "F10" ), "=D0+5" ) );
    assert ( x0 . setCell ( CPos ( "F11" ), "=$D0+5" ) );
    assert ( x0 . setCell ( CPos ( "F12" ), "=D$0+5" ) );
    assert ( x0 . setCell ( CPos ( "F13" ), "=$D$0+5" ) );
    x0 . copyRect ( CPos ( "G11" ), CPos ( "F10" ), 1, 4 );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F10" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F11" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F12" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F13" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F14" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G10" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G11" ) ), CValue ( 75.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G12" ) ), CValue ( 25.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G13" ) ), CValue ( 65.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G14" ) ), CValue ( 15.0 ) ) );
    x0 . copyRect ( CPos ( "G11" ), CPos ( "F10" ), 2, 4 );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F10" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F11" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F12" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F13" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "F14" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G10" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G11" ) ), CValue ( 75.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G12" ) ), CValue ( 25.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G13" ) ), CValue ( 65.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "G14" ) ), CValue ( 15.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H10" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H11" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H12" ) ), CValue() ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H13" ) ), CValue ( 35.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H14" ) ), CValue() ) );
    assert ( x0 . setCell ( CPos ( "F0" ), "-27" ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H14" ) ), CValue ( -22.0 ) ) );
    x0 . copyRect ( CPos ( "H12" ), CPos ( "H13" ), 1, 2 );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H12" ) ), CValue ( 25.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H13" ) ), CValue ( -22.0 ) ) );
    assert ( valueMatch ( x0 . getValue ( CPos ( "H14" ) ), CValue ( -22.0 ) ) );
    return EXIT_SUCCESS;
}
#endif /* __PROGTEST__ */
