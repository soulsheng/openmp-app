
// mfcView.cpp : CmfcView 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "mfc.h"
#endif

#include "mfcDoc.h"
#include "mfcView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CmfcView

IMPLEMENT_DYNCREATE(CmfcView, CView)

BEGIN_MESSAGE_MAP(CmfcView, CView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
END_MESSAGE_MAP()

// CmfcView 构造/析构

CmfcView::CmfcView()
{
	// TODO: 在此处添加构造代码

}

CmfcView::~CmfcView()
{
}

BOOL CmfcView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CView::PreCreateWindow(cs);
}

// CmfcView 绘制

void CmfcView::OnDraw(CDC* pDC)
{
	CmfcDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: 在此处为本机数据添加绘制代码
	RECT rect;
	GetClientRect(&rect);
	pDC->DrawText(_T("Hello World!"), &rect, 0);
}

void CmfcView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CmfcView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CmfcView 诊断

#ifdef _DEBUG
void CmfcView::AssertValid() const
{
	CView::AssertValid();
}

void CmfcView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CmfcDoc* CmfcView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CmfcDoc)));
	return (CmfcDoc*)m_pDocument;
}
#endif //_DEBUG


// CmfcView 消息处理程序
