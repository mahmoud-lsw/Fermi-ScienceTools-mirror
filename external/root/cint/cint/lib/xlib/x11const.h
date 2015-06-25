/* -*- C++ -*- */
/*************************************************************************
 * Copyright(c) 1995~2005  Masaharu Goto (root-cint@cern.ch)
 *
 * For the licensing terms see the file COPYING
 *
 ************************************************************************/
/***********************************************************************
* x11const.h
*
* XLIB define macros exposed to CINT interpreter
***********************************************************************/

#ifndef G__X11CONST_H
#define G__X11CONST_H

/***********************************************************************
* define macro as constant value 
***********************************************************************/
#ifdef NEVER
#define EOF -1
#define NULL 0
#define G__SHAREDLIB 1
#define _SYS_TYPES_INCLUDED 0
#define _SYS_STDSYMS_INCLUDED 0
#define _HPUX_SOURCE 0
#define _INCLUDE__STDC__ 0
#define _INCLUDE_POSIX_SOURCE 0
#define _INCLUDE_POSIX2_SOURCE 0
#define _INCLUDE_XOPEN_SOURCE 0
#define _INCLUDE_AES_SOURCE 0
#define _INCLUDE_HPUX_SOURCE 0
#define _XPG4 0
#define _SIO 0
#define _DEV_T 0
#define _INO_T 0
#define _MODE_T 0
#define _NLINK_T 0
#define _OFF_T 0
#define _PID_T 0
#define _GID_T 0
#define _UID_T 0
#define _TIME_T 0
#define _SIZE_T 0
#define _SSIZE_T 0
#define _SITE_T 0
#define _CNODE_T 0
#define _CLOCK_T 0
#define _KEY_T 0
#define _CADDR_T 0
#define MAXSUSE 65535
#define _AID_T 0
#define UID_NO_CHANGE -1
#define GID_NO_CHANGE -1
#define PGID_NOT_SET -1
#define SID_NOT_SET -1
#define MAXFUPLIM 2048
#define FD_SETSIZE 2048
#define NFDBITS 32
#endif

#define XlibSpecificationRelease 6

#define X_H 0
#define X_PROTOCOL 11
#define X_PROTOCOL_REVISION 0
#define None 0
#define ParentRelative 1
#define CopyFromParent 0
#define PointerWindow 0
#define InputFocus 1
#define PointerRoot 1
#define AnyPropertyType 0
#define AnyKey 0
#define AnyButton 0
#define AllTemporary 0
#define CurrentTime 0
#define NoSymbol 0
#define NoEventMask 0
#define KeyPressMask 1
#define KeyReleaseMask 2
#define ButtonPressMask 4
#define ButtonReleaseMask 8
#define EnterWindowMask 16
#define LeaveWindowMask 32
#define PointerMotionMask 64
#define PointerMotionHintMask 128
#define Button1MotionMask 256
#define Button2MotionMask 512
#define Button3MotionMask 1024
#define Button4MotionMask 2048
#define Button5MotionMask 4096
#define ButtonMotionMask 8192
#define KeymapStateMask 16384
#define ExposureMask 32768
#define VisibilityChangeMask 65536
#define StructureNotifyMask 131072
#define ResizeRedirectMask 262144
#define SubstructureNotifyMask 524288
#define SubstructureRedirectMask 1048576
#define FocusChangeMask 2097152
#define PropertyChangeMask 4194304
#define ColormapChangeMask 8388608
#define OwnerGrabButtonMask 16777216
#define KeyPress 2
#define KeyRelease 3
#define ButtonPress 4
#define ButtonRelease 5
#define MotionNotify 6
#define EnterNotify 7
#define LeaveNotify 8
#define FocusIn 9
#define FocusOut 10
#define KeymapNotify 11
#define Expose 12
#define GraphicsExpose 13
#define NoExpose 14
#define VisibilityNotify 15
#define CreateNotify 16
#define DestroyNotify 17
#define UnmapNotify 18
#define MapNotify 19
#define MapRequest 20
#define ReparentNotify 21
#define ConfigureNotify 22
#define ConfigureRequest 23
#define GravityNotify 24
#define ResizeRequest 25
#define CirculateNotify 26
#define CirculateRequest 27
#define PropertyNotify 28
#define SelectionClear 29
#define SelectionRequest 30
#define SelectionNotify 31
#define ColormapNotify 32
#define ClientMessage 33
#define MappingNotify 34
#define LASTEvent 35
#define ShiftMask 1
#define LockMask 2
#define ControlMask 4
#define Mod1Mask 8
#define Mod2Mask 16
#define Mod3Mask 32
#define Mod4Mask 64
#define Mod5Mask 128
#define ShiftMapIndex 0
#define LockMapIndex 1
#define ControlMapIndex 2
#define Mod1MapIndex 3
#define Mod2MapIndex 4
#define Mod3MapIndex 5
#define Mod4MapIndex 6
#define Mod5MapIndex 7
#define Button1Mask 256
#define Button2Mask 512
#define Button3Mask 1024
#define Button4Mask 2048
#define Button5Mask 4096
#define AnyModifier 32768
#define Button1 1
#define Button2 2
#define Button3 3
#define Button4 4
#define Button5 5
#define NotifyNormal 0
#define NotifyGrab 1
#define NotifyUngrab 2
#define NotifyWhileGrabbed 3
#define NotifyHint 1
#define NotifyAncestor 0
#define NotifyVirtual 1
#define NotifyInferior 2
#define NotifyNonlinear 3
#define NotifyNonlinearVirtual 4
#define NotifyPointer 5
#define NotifyPointerRoot 6
#define NotifyDetailNone 7
#define VisibilityUnobscured 0
#define VisibilityPartiallyObscured 1
#define VisibilityFullyObscured 2
#define PlaceOnTop 0
#define PlaceOnBottom 1
#define FamilyInternet 0
#define FamilyDECnet 1
#define FamilyChaos 2
#define PropertyNewValue 0
#define PropertyDelete 1
#define ColormapUninstalled 0
#define ColormapInstalled 1
#define GrabModeSync 0
#define GrabModeAsync 1
#define GrabSuccess 0
#define AlreadyGrabbed 1
#define GrabInvalidTime 2
#define GrabNotViewable 3
#define GrabFrozen 4
#define AsyncPointer 0
#define SyncPointer 1
#define ReplayPointer 2
#define AsyncKeyboard 3
#define SyncKeyboard 4
#define ReplayKeyboard 5
#define AsyncBoth 6
#define SyncBoth 7
#define RevertToNone 0
#define RevertToPointerRoot 1
#define RevertToParent 2
#define Success 0
#define BadRequest 1
#define BadValue 2
#define BadWindow 3
#define BadPixmap 4
#define BadAtom 5
#define BadCursor 6
#define BadFont 7
#define BadMatch 8
#define BadDrawable 9
#define BadAccess 10
#define BadAlloc 11
#define BadColor 12
#define BadGC 13
#define BadIDChoice 14
#define BadName 15
#define BadLength 16
#define BadImplementation 17
#define FirstExtensionError 128
#define LastExtensionError 255
#define InputOutput 1
#define InputOnly 2
#define CWBackPixmap 1
#define CWBackPixel 2
#define CWBorderPixmap 4
#define CWBorderPixel 8
#define CWBitGravity 16
#define CWWinGravity 32
#define CWBackingStore 64
#define CWBackingPlanes 128
#define CWBackingPixel 256
#define CWOverrideRedirect 512
#define CWSaveUnder 1024
#define CWEventMask 2048
#define CWDontPropagate 4096
#define CWColormap 8192
#define CWCursor 16384
#define CWX 1
#define CWY 2
#define CWWidth 4
#define CWHeight 8
#define CWBorderWidth 16
#define CWSibling 32
#define CWStackMode 64
#define ForgetGravity 0
#define NorthWestGravity 1
#define NorthGravity 2
#define NorthEastGravity 3
#define WestGravity 4
#define CenterGravity 5
#define EastGravity 6
#define SouthWestGravity 7
#define SouthGravity 8
#define SouthEastGravity 9
#define StaticGravity 10
#define UnmapGravity 0
#define NotUseful 0
#define WhenMapped 1
#define Always 2
#define IsUnmapped 0
#define IsUnviewable 1
#define IsViewable 2
#define SetModeInsert 0
#define SetModeDelete 1
#define DestroyAll 0
#define RetainPermanent 1
#define RetainTemporary 2
#define Above 0
#define Below 1
#define TopIf 2
#define BottomIf 3
#define Opposite 4
#define RaiseLowest 0
#define LowerHighest 1
#define PropModeReplace 0
#define PropModePrepend 1
#define PropModeAppend 2
#define GXclear 0
#define GXand 1
#define GXandReverse 2
#define GXcopy 3
#define GXandInverted 4
#define GXnoop 5
#define GXxor 6
#define GXor 7
#define GXnor 8
#define GXequiv 9
#define GXinvert 10
#define GXorReverse 11
#define GXcopyInverted 12
#define GXorInverted 13
#define GXnand 14
#define GXset 15
#define LineSolid 0
#define LineOnOffDash 1
#define LineDoubleDash 2
#define CapNotLast 0
#define CapButt 1
#define CapRound 2
#define CapProjecting 3
#define JoinMiter 0
#define JoinRound 1
#define JoinBevel 2
#define FillSolid 0
#define FillTiled 1
#define FillStippled 2
#define FillOpaqueStippled 3
#define EvenOddRule 0
#define WindingRule 1
#define ClipByChildren 0
#define IncludeInferiors 1
#define Unsorted 0
#define YSorted 1
#define YXSorted 2
#define YXBanded 3
#define CoordModeOrigin 0
#define CoordModePrevious 1
#define Complex 0
#define Nonconvex 1
#define Convex 2
#define ArcChord 0
#define ArcPieSlice 1
#define GCFunction 1
#define GCPlaneMask 2
#define GCForeground 4
#define GCBackground 8
#define GCLineWidth 16
#define GCLineStyle 32
#define GCCapStyle 64
#define GCJoinStyle 128
#define GCFillStyle 256
#define GCFillRule 512
#define GCTile 1024
#define GCStipple 2048
#define GCTileStipXOrigin 4096
#define GCTileStipYOrigin 8192
#define GCFont 16384
#define GCSubwindowMode 32768
#define GCGraphicsExposures 65536
#define GCClipXOrigin 131072
#define GCClipYOrigin 262144
#define GCClipMask 524288
#define GCDashOffset 1048576
#define GCDashList 2097152
#define GCArcMode 4194304
#define GCLastBit 22
#define FontLeftToRight 0
#define FontRightToLeft 1
#define FontChange 255
#define XYBitmap 0
#define XYPixmap 1
#define ZPixmap 2
#define AllocNone 0
#define AllocAll 1
#define DoRed 1
#define DoGreen 2
#define DoBlue 4
#define CursorShape 0
#define TileShape 1
#define StippleShape 2
#define AutoRepeatModeOff 0
#define AutoRepeatModeOn 1
#define AutoRepeatModeDefault 2
#define LedModeOff 0
#define LedModeOn 1
#define KBKeyClickPercent 1
#define KBBellPercent 2
#define KBBellPitch 4
#define KBBellDuration 8
#define KBLed 16
#define KBLedMode 32
#define KBKey 64
#define KBAutoRepeatMode 128
#define MappingSuccess 0
#define MappingBusy 1
#define MappingFailed 2
#define MappingModifier 0
#define MappingKeyboard 1
#define MappingPointer 2
#define DontPreferBlanking 0
#define PreferBlanking 1
#define DefaultBlanking 2
#define DisableScreenSaver 0
#define DisableScreenInterval 0
#define DontAllowExposures 0
#define AllowExposures 1
#define DefaultExposures 2
#define ScreenSaverReset 0
#define ScreenSaverActive 1
#define HostInsert 0
#define HostDelete 1
#define EnableAccess 1
#define DisableAccess 0
#define StaticGray 0
#define GrayScale 1
#define StaticColor 2
#define PseudoColor 3
#define TrueColor 4
#define DirectColor 5
#define LSBFirst 0
#define MSBFirst 1
#define _XFUNCPROTO_H_ 0
#define NeedFunctionPrototypes 0
#define NeedVarargsPrototypes 0
#define _XFUNCPROTOBEGIN 0
#define _XFUNCPROTOEND 0
#define _XOSDEFS_H_ 0
#define G__STDDEF_H 0
#define True 1
#define False 0
#define QueuedAlready 0
#define QueuedAfterReading 1
#define QueuedAfterFlush 2
#define AllPlanes -1
#define XNRequiredCharSet 1074607416
#define XNQueryOrientation 1074607464
#define XNBaseFontName 1074607512
#define XNOMAutomatic 1074607552
#define XNMissingCharSet 1074607592
#define XNDefaultString 1074607632
#define XNOrientation 1074607672
#define XNDirectionalDependentDrawing 1074607712
#define XNContextualDrawing 1074609488
#define XNFontInfo 1074609536
#define XIMPreeditArea 1
#define XIMPreeditCallbacks 2
#define XIMPreeditPosition 4
#define XIMPreeditNothing 8
#define XIMPreeditNone 16
#define XIMStatusArea 256
#define XIMStatusCallbacks 512
#define XIMStatusNothing 1024
#define XIMStatusNone 2048
#define XNVaNestedList 1074658872
#define XNQueryInputStyle 1074658912
#define XNClientWindow 1074658960
#define XNInputStyle 1074659000
#define XNFocusWindow 1074660760
#define XNResourceName 1074660800
#define XNResourceClass 1074660840
#define XNGeometryCallback 1074660880
#define XNDestroyCallback 1074660928
#define XNFilterEvents 1074660976
#define XNPreeditStartCallback 1074661016
#define XNPreeditDoneCallback 1074661064
#define XNPreeditDrawCallback 1074661112
#define XNPreeditCaretCallback 1074661160
#define XNPreeditAttributes 1074662928
#define XNStatusStartCallback 1074662976
#define XNStatusDoneCallback 1074663024
#define XNStatusDrawCallback 1074663072
#define XNStatusAttributes 1074663120
#define XNArea 1074663168
#define XNAreaNeeded 1074663200
#define XNSpotLocation 1074663240
#define XNColormap 1074663280
#define XNStdColormap 1074663320
#define XNForeground 1074665080
#define XNBackground 1074665120
#define XNBackgroundPixmap 1074665160
#define XNFontSet 1074665208
#define XNLineSpace 1074665248
#define XNCursor 1074665288
#define XNQueryIMValuesList 1074665320
#define XNQueryICValuesList 1074665368
#define XNVisiblePosition 1074665416
#define XNR6PreeditCallback 1074665464
#define XNStringConversionCallback 1074667232
#define XNStringConversion 1074667288
#define XNResetState 1074667336
#define XNHotKey 1074667376
#define XNHotKeyState 1074667408
#define XNPreeditState 1074667448
#define XNSeparatorofNestedList 1074667488
#define XBufferOverflow -1
#define XLookupNone 1
#define XLookupChars 2
#define XLookupKeySym 3
#define XLookupBoth 4
#define XIMReverse 1
#define XIMUnderline 2
#define XIMHighlight 4
#define XIMPrimary 32
#define XIMSecondary 64
#define XIMTertiary 128
#define XIMVisibleToForward 256
#define XIMVisibleToBackword 512
#define XIMVisibleToCenter 1024
#define XIMPreeditUnKnown 0
#define XIMPreeditEnable 1
#define XIMPreeditDisable 2
#define XIMInitialState 1
#define XIMPreserveState 2
#define XIMStringConversionLeftEdge 1
#define XIMStringConversionRightEdge 2
#define XIMStringConversionTopEdge 4
#define XIMStringConversionBottomEdge 8
#define XIMStringConversionConcealed 16
#define XIMStringConversionWrapped 32
#define XIMStringConversionBuffer 1
#define XIMStringConversionLine 2
#define XIMStringConversionWord 3
#define XIMStringConversionChar 4
#define XIMStringConversionSubstitution 1
#define XIMStringConversionRetrival 2
#define XIMHotKeyStateON 1
#define XIMHotKeyStateOFF 2
#define _XUTIL_H_ 0
#define NoValue 0
#define XValue 1
#define YValue 2
#define WidthValue 4
#define HeightValue 8
#define AllValues 15
#define XNegative 16
#define YNegative 32
#define USPosition 1
#define USSize 2
#define PPosition 4
#define PSize 8
#define PMinSize 16
#define PMaxSize 32
#define PResizeInc 64
#define PAspect 128
#define PBaseSize 256
#define PWinGravity 512
#define PAllHints 252
#define InputHint 1
#define StateHint 2
#define IconPixmapHint 4
#define IconWindowHint 8
#define IconPositionHint 16
#define IconMaskHint 32
#define WindowGroupHint 64
#define AllHints 111
#define XUrgencyHint 256
#define WithdrawnState 0
#define NormalState 1
#define IconicState 3
#define DontCareState 0
#define ZoomState 2
#define InactiveState 4
#define XNoMemory -1
#define XLocaleNotSupported -2
#define XConverterNotFound -3
#define RectangleOut 0
#define RectangleIn 1
#define RectanglePart 2
#define VisualNoMask 0
#define VisualIDMask 1
#define VisualScreenMask 2
#define VisualDepthMask 4
#define VisualClassMask 8
#define VisualRedMaskMask 16
#define VisualGreenMaskMask 32
#define VisualBlueMaskMask 64
#define VisualColormapSizeMask 128
#define VisualBitsPerRGBMask 256
#define VisualAllMask 511
#define ReleaseByFreeingColormap 1
#define BitmapSuccess 0
#define BitmapOpenFailed 1
#define BitmapFileInvalid 2
#define BitmapNoMemory 3
#define XCSUCCESS 0
#define XCNOMEM 1
#define XCNOENT 2
#define _XOS_H_ 0
#define index 1073813410
#define rindex 1073813818
#define G__FCNTL_H 0
#define _SYS_TIME_INCLUDED 0
#define CLOCKS_PER_SEC 1000000
#define _STRUCT_TIMEVAL 0
#define DST_NONE 0
#define DST_USA 1
#define DST_AUST 2
#define DST_WET 3
#define DST_MET 4
#define DST_EET 5
#define ITIMER_REAL 0
#define ITIMER_VIRTUAL 1
#define ITIMER_PROF 2
#define XK_MISCELLANY 0
#define XK_XKB_KEYS 0
#define XK_LATIN1 0
#define XK_LATIN2 0
#define XK_LATIN3 0
#define XK_LATIN4 0
#define XK_GREEK 0
#define XK_VoidSymbol 16777215
#define XK_BackSpace 65288
#define XK_Tab 65289
#define XK_Linefeed 65290
#define XK_Clear 65291
#define XK_Return 65293
#define XK_Pause 65299
#define XK_Scroll_Lock 65300
#define XK_Sys_Req 65301
#define XK_Escape 65307
#define XK_Delete 65535
#define XK_Multi_key 65312
#define XK_Kanji 65313
#define XK_Muhenkan 65314
#define XK_Henkan_Mode 65315
#define XK_Henkan 65315
#define XK_Romaji 65316
#define XK_Hiragana 65317
#define XK_Katakana 65318
#define XK_Hiragana_Katakana 65319
#define XK_Zenkaku 65320
#define XK_Hankaku 65321
#define XK_Zenkaku_Hankaku 65322
#define XK_Touroku 65323
#define XK_Massyo 65324
#define XK_Kana_Lock 65325
#define XK_Kana_Shift 65326
#define XK_Eisu_Shift 65327
#define XK_Eisu_toggle 65328
#define XK_Home 65360
#define XK_Left 65361
#define XK_Up 65362
#define XK_Right 65363
#define XK_Down 65364
#define XK_Prior 65365
#define XK_Page_Up 65365
#define XK_Next 65366
#define XK_Page_Down 65366
#define XK_End 65367
#define XK_Begin 65368
#define XK_Select 65376
#define XK_Print 65377
#define XK_Execute 65378
#define XK_Insert 65379
#define XK_Undo 65381
#define XK_Redo 65382
#define XK_Menu 65383
#define XK_Find 65384
#define XK_Cancel 65385
#define XK_Help 65386
#define XK_Break 65387
#define XK_Mode_switch 65406
#define XK_script_switch 65406
#define XK_Num_Lock 65407
#define XK_KP_Space 65408
#define XK_KP_Tab 65417
#define XK_KP_Enter 65421
#define XK_KP_F1 65425
#define XK_KP_F2 65426
#define XK_KP_F3 65427
#define XK_KP_F4 65428
#define XK_KP_Home 65429
#define XK_KP_Left 65430
#define XK_KP_Up 65431
#define XK_KP_Right 65432
#define XK_KP_Down 65433
#define XK_KP_Prior 65434
#define XK_KP_Page_Up 65434
#define XK_KP_Next 65435
#define XK_KP_Page_Down 65435
#define XK_KP_End 65436
#define XK_KP_Begin 65437
#define XK_KP_Insert 65438
#define XK_KP_Delete 65439
#define XK_KP_Equal 65469
#define XK_KP_Multiply 65450
#define XK_KP_Add 65451
#define XK_KP_Separator 65452
#define XK_KP_Subtract 65453
#define XK_KP_Decimal 65454
#define XK_KP_Divide 65455
#define XK_KP_0 65456
#define XK_KP_1 65457
#define XK_KP_2 65458
#define XK_KP_3 65459
#define XK_KP_4 65460
#define XK_KP_5 65461
#define XK_KP_6 65462
#define XK_KP_7 65463
#define XK_KP_8 65464
#define XK_KP_9 65465
#define XK_F1 65470
#define XK_F2 65471
#define XK_F3 65472
#define XK_F4 65473
#define XK_F5 65474
#define XK_F6 65475
#define XK_F7 65476
#define XK_F8 65477
#define XK_F9 65478
#define XK_F10 65479
#define XK_F11 65480
#define XK_L1 65480
#define XK_F12 65481
#define XK_L2 65481
#define XK_F13 65482
#define XK_L3 65482
#define XK_F14 65483
#define XK_L4 65483
#define XK_F15 65484
#define XK_L5 65484
#define XK_F16 65485
#define XK_L6 65485
#define XK_F17 65486
#define XK_L7 65486
#define XK_F18 65487
#define XK_L8 65487
#define XK_F19 65488
#define XK_L9 65488
#define XK_F20 65489
#define XK_L10 65489
#define XK_F21 65490
#define XK_R1 65490
#define XK_F22 65491
#define XK_R2 65491
#define XK_F23 65492
#define XK_R3 65492
#define XK_F24 65493
#define XK_R4 65493
#define XK_F25 65494
#define XK_R5 65494
#define XK_F26 65495
#define XK_R6 65495
#define XK_F27 65496
#define XK_R7 65496
#define XK_F28 65497
#define XK_R8 65497
#define XK_F29 65498
#define XK_R9 65498
#define XK_F30 65499
#define XK_R10 65499
#define XK_F31 65500
#define XK_R11 65500
#define XK_F32 65501
#define XK_R12 65501
#define XK_F33 65502
#define XK_R13 65502
#define XK_F34 65503
#define XK_R14 65503
#define XK_F35 65504
#define XK_R15 65504
#define XK_Shift_L 65505
#define XK_Shift_R 65506
#define XK_Control_L 65507
#define XK_Control_R 65508
#define XK_Caps_Lock 65509
#define XK_Shift_Lock 65510
#define XK_Meta_L 65511
#define XK_Meta_R 65512
#define XK_Alt_L 65513
#define XK_Alt_R 65514
#define XK_Super_L 65515
#define XK_Super_R 65516
#define XK_Hyper_L 65517
#define XK_Hyper_R 65518
#define XK_ISO_Lock 65025
#define XK_ISO_Level2_Latch 65026
#define XK_ISO_Level3_Shift 65027
#define XK_ISO_Level3_Latch 65028
#define XK_ISO_Level3_Lock 65029
#define XK_ISO_Group_Shift 65406
#define XK_ISO_Group_Latch 65030
#define XK_ISO_Group_Lock 65031
#define XK_ISO_Next_Group 65032
#define XK_ISO_Next_Group_Lock 65033
#define XK_ISO_Prev_Group 65034
#define XK_ISO_Prev_Group_Lock 65035
#define XK_ISO_First_Group 65036
#define XK_ISO_First_Group_Lock 65037
#define XK_ISO_Last_Group 65038
#define XK_ISO_Last_Group_Lock 65039
#define XK_ISO_Left_Tab 65056
#define XK_ISO_Move_Line_Up 65057
#define XK_ISO_Move_Line_Down 65058
#define XK_ISO_Partial_Line_Up 65059
#define XK_ISO_Partial_Line_Down 65060
#define XK_ISO_Partial_Space_Left 65061
#define XK_ISO_Partial_Space_Right 65062
#define XK_ISO_Set_Margin_Left 65063
#define XK_ISO_Set_Margin_Right 65064
#define XK_ISO_Release_Margin_Left 65065
#define XK_ISO_Release_Margin_Right 65066
#define XK_ISO_Release_Both_Margins 65067
#define XK_ISO_Fast_Cursor_Left 65068
#define XK_ISO_Fast_Cursor_Right 65069
#define XK_ISO_Fast_Cursor_Up 65070
#define XK_ISO_Fast_Cursor_Down 65071
#define XK_ISO_Continuous_Underline 65072
#define XK_ISO_Discontinuous_Underline 65073
#define XK_ISO_Emphasize 65074
#define XK_ISO_Center_Object 65075
#define XK_ISO_Enter 65076
#define XK_dead_grave 65104
#define XK_dead_acute 65105
#define XK_dead_circumflex 65106
#define XK_dead_tilde 65107
#define XK_dead_macron 65108
#define XK_dead_breve 65109
#define XK_dead_abovedot 65110
#define XK_dead_diaeresis 65111
#define XK_dead_abovering 65112
#define XK_dead_doubleacute 65113
#define XK_dead_caron 65114
#define XK_dead_cedilla 65115
#define XK_dead_ogonek 65116
#define XK_dead_iota 65117
#define XK_dead_voiced_sound 65118
#define XK_dead_semivoiced_sound 65119
#define XK_First_Virtual_Screen 65232
#define XK_Prev_Virtual_Screen 65233
#define XK_Next_Virtual_Screen 65234
#define XK_Last_Virtual_Screen 65236
#define XK_Terminate_Server 65237
#define XK_Pointer_Left 65248
#define XK_Pointer_Right 65249
#define XK_Pointer_Up 65250
#define XK_Pointer_Down 65251
#define XK_Pointer_UpLeft 65252
#define XK_Pointer_UpRight 65253
#define XK_Pointer_DownLeft 65254
#define XK_Pointer_DownRight 65255
#define XK_Pointer_Button_Dflt 65256
#define XK_Pointer_Button1 65257
#define XK_Pointer_Button2 65258
#define XK_Pointer_Button3 65259
#define XK_Pointer_Button4 65260
#define XK_Pointer_Button5 65261
#define XK_Pointer_DblClick_Dflt 65262
#define XK_Pointer_DblClick1 65263
#define XK_Pointer_DblClick2 65264
#define XK_Pointer_DblClick3 65265
#define XK_Pointer_DblClick4 65266
#define XK_Pointer_DblClick5 65267
#define XK_Pointer_Drag_Dflt 65268
#define XK_Pointer_Drag1 65269
#define XK_Pointer_Drag2 65270
#define XK_Pointer_Drag3 65271
#define XK_Pointer_Drag4 65272
#define XK_Pointer_EnableKeys 65273
#define XK_Pointer_Accelerate 65274
#define XK_Pointer_DfltBtnNext 65275
#define XK_Pointer_DfltBtnPrev 65276
#define XK_space 32
#define XK_exclam 33
#define XK_quotedbl 34
#define XK_numbersign 35
#define XK_dollar 36
#define XK_percent 37
#define XK_ampersand 38
#define XK_apostrophe 39
#define XK_quoteright 39
#define XK_parenleft 40
#define XK_parenright 41
#define XK_asterisk 42
#define XK_plus 43
#define XK_comma 44
#define XK_minus 45
#define XK_period 46
#define XK_slash 47
#define XK_0 48
#define XK_1 49
#define XK_2 50
#define XK_3 51
#define XK_4 52
#define XK_5 53
#define XK_6 54
#define XK_7 55
#define XK_8 56
#define XK_9 57
#define XK_colon 58
#define XK_semicolon 59
#define XK_less 60
#define XK_equal 61
#define XK_greater 62
#define XK_question 63
#define XK_at 64
#define XK_A 65
#define XK_B 66
#define XK_C 67
#define XK_D 68
#define XK_E 69
#define XK_F 70
#define XK_G 71
#define XK_H 72
#define XK_I 73
#define XK_J 74
#define XK_K 75
#define XK_L 76
#define XK_M 77
#define XK_N 78
#define XK_O 79
#define XK_P 80
#define XK_Q 81
#define XK_R 82
#define XK_S 83
#define XK_T 84
#define XK_U 85
#define XK_V 86
#define XK_W 87
#define XK_X 88
#define XK_Y 89
#define XK_Z 90
#define XK_bracketleft 91
#define XK_backslash 92
#define XK_bracketright 93
#define XK_asciicircum 94
#define XK_underscore 95
#define XK_grave 96
#define XK_quoteleft 96
#define XK_a 97
#define XK_b 98
#define XK_c 99
#define XK_d 100
#define XK_e 101
#define XK_f 102
#define XK_g 103
#define XK_h 104
#define XK_i 105
#define XK_j 106
#define XK_k 107
#define XK_l 108
#define XK_m 109
#define XK_n 110
#define XK_o 111
#define XK_p 112
#define XK_q 113
#define XK_r 114
#define XK_s 115
#define XK_t 116
#define XK_u 117
#define XK_v 118
#define XK_w 119
#define XK_x 120
#define XK_y 121
#define XK_z 122
#define XK_braceleft 123
#define XK_bar 124
#define XK_braceright 125
#define XK_asciitilde 126
#define XK_nobreakspace 160
#define XK_exclamdown 161
#define XK_cent 162
#define XK_sterling 163
#define XK_currency 164
#define XK_yen 165
#define XK_brokenbar 166
#define XK_section 167
#define XK_diaeresis 168
#define XK_copyright 169
#define XK_ordfeminine 170
#define XK_guillemotleft 171
#define XK_notsign 172
#define XK_hyphen 173
#define XK_registered 174
#define XK_macron 175
#define XK_degree 176
#define XK_plusminus 177
#define XK_twosuperior 178
#define XK_threesuperior 179
#define XK_acute 180
#define XK_mu 181
#define XK_paragraph 182
#define XK_periodcentered 183
#define XK_cedilla 184
#define XK_onesuperior 185
#define XK_masculine 186
#define XK_guillemotright 187
#define XK_onequarter 188
#define XK_onehalf 189
#define XK_threequarters 190
#define XK_questiondown 191
#define XK_Agrave 192
#define XK_Aacute 193
#define XK_Acircumflex 194
#define XK_Atilde 195
#define XK_Adiaeresis 196
#define XK_Aring 197
#define XK_AE 198
#define XK_Ccedilla 199
#define XK_Egrave 200
#define XK_Eacute 201
#define XK_Ecircumflex 202
#define XK_Ediaeresis 203
#define XK_Igrave 204
#define XK_Iacute 205
#define XK_Icircumflex 206
#define XK_Idiaeresis 207
#define XK_ETH 208
#define XK_Eth 208
#define XK_Ntilde 209
#define XK_Ograve 210
#define XK_Oacute 211
#define XK_Ocircumflex 212
#define XK_Otilde 213
#define XK_Odiaeresis 214
#define XK_multiply 215
#define XK_Ooblique 216
#define XK_Ugrave 217
#define XK_Uacute 218
#define XK_Ucircumflex 219
#define XK_Udiaeresis 220
#define XK_Yacute 221
#define XK_THORN 222
#define XK_Thorn 222
#define XK_ssharp 223
#define XK_agrave 224
#define XK_aacute 225
#define XK_acircumflex 226
#define XK_atilde 227
#define XK_adiaeresis 228
#define XK_aring 229
#define XK_ae 230
#define XK_ccedilla 231
#define XK_egrave 232
#define XK_eacute 233
#define XK_ecircumflex 234
#define XK_ediaeresis 235
#define XK_igrave 236
#define XK_iacute 237
#define XK_icircumflex 238
#define XK_idiaeresis 239
#define XK_eth 240
#define XK_ntilde 241
#define XK_ograve 242
#define XK_oacute 243
#define XK_ocircumflex 244
#define XK_otilde 245
#define XK_odiaeresis 246
#define XK_division 247
#define XK_oslash 248
#define XK_ugrave 249
#define XK_uacute 250
#define XK_ucircumflex 251
#define XK_udiaeresis 252
#define XK_yacute 253
#define XK_thorn 254
#define XK_ydiaeresis 255
#define XK_Aogonek 417
#define XK_breve 418
#define XK_Lstroke 419
#define XK_Lcaron 421
#define XK_Sacute 422
#define XK_Scaron 425
#define XK_Scedilla 426
#define XK_Tcaron 427
#define XK_Zacute 428
#define XK_Zcaron 430
#define XK_Zabovedot 431
#define XK_aogonek 433
#define XK_ogonek 434
#define XK_lstroke 435
#define XK_lcaron 437
#define XK_sacute 438
#define XK_caron 439
#define XK_scaron 441
#define XK_scedilla 442
#define XK_tcaron 443
#define XK_zacute 444
#define XK_doubleacute 445
#define XK_zcaron 446
#define XK_zabovedot 447
#define XK_Racute 448
#define XK_Abreve 451
#define XK_Lacute 453
#define XK_Cacute 454
#define XK_Ccaron 456
#define XK_Eogonek 458
#define XK_Ecaron 460
#define XK_Dcaron 463
#define XK_Dstroke 464
#define XK_Nacute 465
#define XK_Ncaron 466
#define XK_Odoubleacute 469
#define XK_Rcaron 472
#define XK_Uring 473
#define XK_Udoubleacute 475
#define XK_Tcedilla 478
#define XK_racute 480
#define XK_abreve 483
#define XK_lacute 485
#define XK_cacute 486
#define XK_ccaron 488
#define XK_eogonek 490
#define XK_ecaron 492
#define XK_dcaron 495
#define XK_dstroke 496
#define XK_nacute 497
#define XK_ncaron 498
#define XK_odoubleacute 501
#define XK_udoubleacute 507
#define XK_rcaron 504
#define XK_uring 505
#define XK_tcedilla 510
#define XK_abovedot 511
#define XK_Hstroke 673
#define XK_Hcircumflex 678
#define XK_Iabovedot 681
#define XK_Gbreve 683
#define XK_Jcircumflex 684
#define XK_hstroke 689
#define XK_hcircumflex 694
#define XK_idotless 697
#define XK_gbreve 699
#define XK_jcircumflex 700
#define XK_Cabovedot 709
#define XK_Ccircumflex 710
#define XK_Gabovedot 725
#define XK_Gcircumflex 728
#define XK_Ubreve 733
#define XK_Scircumflex 734
#define XK_cabovedot 741
#define XK_ccircumflex 742
#define XK_gabovedot 757
#define XK_gcircumflex 760
#define XK_ubreve 765
#define XK_scircumflex 766
#define XK_kra 930
#define XK_kappa 930
#define XK_Rcedilla 931
#define XK_Itilde 933
#define XK_Lcedilla 934
#define XK_Emacron 938
#define XK_Gcedilla 939
#define XK_Tslash 940
#define XK_rcedilla 947
#define XK_itilde 949
#define XK_lcedilla 950
#define XK_emacron 954
#define XK_gcedilla 955
#define XK_tslash 956
#define XK_ENG 957
#define XK_eng 959
#define XK_Amacron 960
#define XK_Iogonek 967
#define XK_Eabovedot 972
#define XK_Imacron 975
#define XK_Ncedilla 977
#define XK_Omacron 978
#define XK_Kcedilla 979
#define XK_Uogonek 985
#define XK_Utilde 989
#define XK_Umacron 990
#define XK_amacron 992
#define XK_iogonek 999
#define XK_eabovedot 1004
#define XK_imacron 1007
#define XK_ncedilla 1009
#define XK_omacron 1010
#define XK_kcedilla 1011
#define XK_uogonek 1017
#define XK_utilde 1021
#define XK_umacron 1022
#define XK_Greek_ALPHAaccent 1953
#define XK_Greek_EPSILONaccent 1954
#define XK_Greek_ETAaccent 1955
#define XK_Greek_IOTAaccent 1956
#define XK_Greek_IOTAdiaeresis 1957
#define XK_Greek_OMICRONaccent 1959
#define XK_Greek_UPSILONaccent 1960
#define XK_Greek_UPSILONdieresis 1961
#define XK_Greek_OMEGAaccent 1963
#define XK_Greek_accentdieresis 1966
#define XK_Greek_horizbar 1967
#define XK_Greek_alphaaccent 1969
#define XK_Greek_epsilonaccent 1970
#define XK_Greek_etaaccent 1971
#define XK_Greek_iotaaccent 1972
#define XK_Greek_iotadieresis 1973
#define XK_Greek_iotaaccentdieresis 1974
#define XK_Greek_omicronaccent 1975
#define XK_Greek_upsilonaccent 1976
#define XK_Greek_upsilondieresis 1977
#define XK_Greek_upsilonaccentdieresis 1978
#define XK_Greek_omegaaccent 1979
#define XK_Greek_ALPHA 1985
#define XK_Greek_BETA 1986
#define XK_Greek_GAMMA 1987
#define XK_Greek_DELTA 1988
#define XK_Greek_EPSILON 1989
#define XK_Greek_ZETA 1990
#define XK_Greek_ETA 1991
#define XK_Greek_THETA 1992
#define XK_Greek_IOTA 1993
#define XK_Greek_KAPPA 1994
#define XK_Greek_LAMDA 1995
#define XK_Greek_LAMBDA 1995
#define XK_Greek_MU 1996
#define XK_Greek_NU 1997
#define XK_Greek_XI 1998
#define XK_Greek_OMICRON 1999
#define XK_Greek_PI 2000
#define XK_Greek_RHO 2001
#define XK_Greek_SIGMA 2002
#define XK_Greek_TAU 2004
#define XK_Greek_UPSILON 2005
#define XK_Greek_PHI 2006
#define XK_Greek_CHI 2007
#define XK_Greek_PSI 2008
#define XK_Greek_OMEGA 2009
#define XK_Greek_alpha 2017
#define XK_Greek_beta 2018
#define XK_Greek_gamma 2019
#define XK_Greek_delta 2020
#define XK_Greek_epsilon 2021
#define XK_Greek_zeta 2022
#define XK_Greek_eta 2023
#define XK_Greek_theta 2024
#define XK_Greek_iota 2025
#define XK_Greek_kappa 2026
#define XK_Greek_lamda 2027
#define XK_Greek_lambda 2027
#define XK_Greek_mu 2028
#define XK_Greek_nu 2029
#define XK_Greek_xi 2030
#define XK_Greek_omicron 2031
#define XK_Greek_pi 2032
#define XK_Greek_rho 2033
#define XK_Greek_sigma 2034
#define XK_Greek_finalsmallsigma 2035
#define XK_Greek_tau 2036
#define XK_Greek_upsilon 2037
#define XK_Greek_phi 2038
#define XK_Greek_chi 2039
#define XK_Greek_psi 2040
#define XK_Greek_omega 2041
#define XK_Greek_switch 65406

#endif 

